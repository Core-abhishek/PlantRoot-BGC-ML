#!/usr/bin/env python3
import csv
import io
import json
import re
import tarfile
import zipfile
from collections import Counter
from pathlib import Path

TAR_PATH = Path("raw/zenodo/CRBC_BGC_part_1.tar.gz")
JSON_DIR = Path("processed/bgc_jsons")
OUT_CSV = Path("results/03_bgc_extraction/all_bgc_summary.csv")
LIMIT = 200

def parse_location(loc_str):
    match = re.search(r"\[(\d+):(\d+)\]", str(loc_str))
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None

def extract_jsons():
    JSON_DIR.mkdir(parents=True, exist_ok=True)
    extracted = 0
    with tarfile.open(TAR_PATH, "r:gz") as tar:
        for member in tar.getmembers():
            if extracted >= LIMIT:
                break
            if not member.isfile() or not member.name.endswith(".zip"):
                continue

            genome_id = Path(member.name).stem
            out_json = JSON_DIR / f"{genome_id}.json"
            if out_json.exists():
                print(f"exists {genome_id}")
                extracted += 1
                continue

            zf = tar.extractfile(member)
            if zf is None:
                print(f"skip {member.name}")
                continue

            with zipfile.ZipFile(io.BytesIO(zf.read())) as z:
                json_names = [n for n in z.namelist() if n.endswith(".json")]
                if not json_names:
                    print(f"no json {member.name}")
                    continue
                out_json.write_bytes(z.read(json_names[0]))

            extracted += 1
            if extracted % 10 == 0:
                print(f"extracted {extracted}")

    print(f"json extraction done: {extracted}")

def summarize_jsons():
    rows = []
    json_files = sorted(JSON_DIR.glob("*.json"))
    print(f"summarizing {len(json_files)} JSON files")

    for json_file in json_files:
        genome_id = json_file.stem
        try:
            with open(json_file) as f:
                data = json.load(f)
        except Exception as e:
            print(f"error reading {json_file.name}: {e}")
            continue

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

    with open(OUT_CSV, "w", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile,
            fieldnames=[
                "genome_id", "record_id", "products", "start", "end",
                "core_start", "core_end", "knowncluster_hits",
                "func_transport", "func_regulatory", "func_biosynthetic",
                "func_biosynthetic_additional", "func_other",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"extracted {len(rows)} BGC regions from {len(json_files)} genomes")
    print(f"saved to: {OUT_CSV}")

if __name__ == "__main__":
    extract_jsons()
    summarize_jsons()
