#!/usr/bin/env python3
import csv
import io
import json
import re
import tarfile
import zipfile
from collections import Counter
from pathlib import Path

TAR_FILES = [
    Path("raw/zenodo/CRBC_BGC_part_1.tar.gz"),
    Path("raw/zenodo/CRBC_BGC_part_2.tar.gz"),
]
JSON_DIR = Path("processed/bgc_jsons")
OUT_CSV = Path("results/03_bgc_extraction/all_bgc_summary.csv")
DONE_FILE = Path("results/03_bgc_extraction/processed_genomes.txt")
BATCH_SIZE = 200

FIELDNAMES = [
    "genome_id",
    "record_id",
    "products",
    "start",
    "end",
    "core_start",
    "core_end",
    "knowncluster_hits",
    "func_transport",
    "func_regulatory",
    "func_biosynthetic",
    "func_biosynthetic_additional",
    "func_other",
]

def parse_location(loc_str):
    match = re.search(r"\[(\d+):(\d+)\]", str(loc_str))
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None

def load_done():
    done = set()
    if DONE_FILE.exists():
        done.update(DONE_FILE.read_text().split())
    if OUT_CSV.exists():
        with open(OUT_CSV, newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames and "genome_id" in reader.fieldnames:
                for row in reader:
                    if row.get("genome_id"):
                        done.add(row["genome_id"])
    return done

def append_done(genome_ids):
    DONE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(DONE_FILE, "a") as f:
        for genome_id in genome_ids:
            f.write(genome_id + "\n")

def write_rows(rows, write_header):
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if write_header else "a"
    with open(OUT_CSV, mode, newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)

def rows_from_json(json_path):
    genome_id = json_path.stem
    rows = []
    with open(json_path) as f:
        data = json.load(f)

    for record in data.get("records", []):
        record_id = record.get("id", "unknown")
        modules = record.get("modules", {})

        gene_to_function = {}
        genefunc = modules.get("antismash.detection.genefunctions", {})
        for tool in genefunc.get("tools", []):
            if tool.get("tool") == "smcogs":
                gene_to_function = tool.get("mapping", {})

        gene_locations = {}
        for feat in record.get("features", []):
            if feat.get("type") != "CDS":
                continue
            quals = feat.get("qualifiers", {})
            gene_id = (quals.get("gene") or quals.get("ID") or [None])[0]
            if not gene_id:
                continue
            g_start, g_end = parse_location(feat.get("location"))
            if g_start is not None:
                gene_locations[gene_id] = (g_start, g_end)

        knowncluster_hits = 0
        known = modules.get("antismash.modules.clusterblast", {}).get("knowncluster", {})
        results_list = known.get("results", [])
        if results_list:
            knowncluster_hits = results_list[0].get("total_hits", 0)

        for area in record.get("areas", []):
            area_start = area.get("start")
            area_end = area.get("end")
            products = area.get("products", [])
            protoclusters = area.get("protoclusters", {})
            core_start = core_end = product_from_proto = None
            if protoclusters:
                first_proto = list(protoclusters.values())[0]
                core_start = first_proto.get("core_start")
                core_end = first_proto.get("core_end")
                product_from_proto = first_proto.get("product")

            function_counts = Counter()
            if area_start is not None and area_end is not None:
                for gene_id, (g_start, g_end) in gene_locations.items():
                    if g_start < area_end and g_end > area_start:
                        func = gene_to_function.get(gene_id)
                        if func:
                            function_counts[func] += 1

            rows.append({
                "genome_id": genome_id,
                "record_id": record_id,
                "products": ";".join(products) if products else (product_from_proto or ""),
                "start": area_start,
                "end": area_end,
                "core_start": core_start,
                "core_end": core_end,
                "knowncluster_hits": knowncluster_hits,
                "func_transport": function_counts.get("transport", 0),
                "func_regulatory": function_counts.get("regulatory", 0),
                "func_biosynthetic": function_counts.get("biosynthetic", 0),
                "func_biosynthetic_additional": function_counts.get("biosynthetic-additional", 0),
                "func_other": function_counts.get("other", 0),
            })
    return rows

def extract_one_json(tar, member, out_json):
    zf = tar.extractfile(member)
    if zf is None:
        return False
    with zipfile.ZipFile(io.BytesIO(zf.read())) as z:
        json_names = [n for n in z.namelist() if n.endswith(".json")]
        if not json_names:
            return False
        out_json.write_bytes(z.read(json_names[0]))
    return True

def main():
    JSON_DIR.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    done = load_done()
    write_header = not OUT_CSV.exists()
    total_genomes = 0
    total_rows = 0

    for tar_path in TAR_FILES:
        if not tar_path.exists():
            print(f"missing {tar_path}")
            continue

        print(f"opening {tar_path}")
        with tarfile.open(tar_path, "r:gz") as tar:
            members = [
                m for m in tar.getmembers()
                if m.isfile() and m.name.endswith(".zip")
            ]
            i = 0
            while i < len(members):
                batch_ids = []
                while i < len(members) and len(batch_ids) < BATCH_SIZE:
                    member = members[i]
                    i += 1
                    genome_id = Path(member.name).stem
                    if genome_id in done:
                        continue
                    out_json = JSON_DIR / f"{genome_id}.json"
                    ok = extract_one_json(tar, member, out_json)
                    if not ok:
                        print(f"no json {member.name}")
                        append_done([genome_id])
                        done.add(genome_id)
                        continue
                    batch_ids.append(genome_id)

                if not batch_ids:
                    continue

                batch_rows = []
                for genome_id in batch_ids:
                    json_path = JSON_DIR / f"{genome_id}.json"
                    try:
                        batch_rows.extend(rows_from_json(json_path))
                    except Exception as e:
                        print(f"error summarizing {genome_id}: {e}")
                    finally:
                        if json_path.exists():
                            json_path.unlink()

                write_rows(batch_rows, write_header)
                write_header = False
                append_done(batch_ids)
                done.update(batch_ids)
                total_genomes += len(batch_ids)
                total_rows += len(batch_rows)
                print(
                    f"batch done: {len(batch_ids)} genomes, "
                    f"{len(batch_rows)} BGC rows, "
                    f"total genomes {len(done)}"
                )

    print(f"finished. new genomes this run: {total_genomes}")
    print(f"new BGC rows this run: {total_rows}")
    print(f"csv: {OUT_CSV}")

if __name__ == "__main__":
    main()
