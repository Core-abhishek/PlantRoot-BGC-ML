#!/usr/bin/env python3
import json
import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

# ========== SETTINGS ==========
input_dir = Path("processed/bgc_jsons")
output_file = Path("results/03_bgc_extraction/all_bgc_summary.csv")
# ==============================

output_file.parent.mkdir(parents=True, exist_ok=True)

def parse_location(loc_str):
    """Parse location string like '[1253039:1275388](-)' → (start, end)"""
    if not loc_str:
        return None, None
    match = re.search(r'\[(\d+):(\d+)\]', str(loc_str))
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None

results = []
json_files = sorted(input_dir.glob("*.json"))
print(f"Found {len(json_files)} JSON files")

for json_file in json_files:
    genome_id = json_file.stem
    
    try:
        with open(json_file) as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_file.name}: {e}")
        continue
    
    for record in data.get("records", []):
        record_id = record.get("id", "unknown")
        modules = record.get("modules", {})
        
        # ----- Build gene → function map (SMCOGs) -----
        gene_to_function = {}
        try:
            genefunc = modules.get("antismash.detection.genefunctions", {})
            for tool in genefunc.get("tools", []):
                if tool.get("tool") == "smcogs":
                    gene_to_function = tool.get("mapping", {})
        except:
            pass
        
        # ----- Build gene → location map from features -----
        gene_locations = {}  # gene_id → (start, end)
        for feat in record.get("features", []):
            if feat.get("type") != "CDS":
                continue
            quals = feat.get("qualifiers", {})
            gene_id = None
            if quals.get("gene"):
                gene_id = quals["gene"][0]
            elif quals.get("ID"):
                gene_id = quals["ID"][0]
            
            if gene_id:
                start, end = parse_location(feat.get("location"))
                if start is not None:
                    gene_locations[gene_id] = (start, end)
        
        # ----- KnownClusterBlast (still contig-level for now) -----
        knowncluster_hits = 0
        try:
            known = modules.get("antismash.modules.clusterblast", {}).get("knowncluster", {})
            results_list = known.get("results", [])
            if results_list:
                knowncluster_hits = results_list[0].get("total_hits", 0)
        except:
            pass
        
        # ----- Process each BGC region -----
        for area in record.get("areas", []):
            products = area.get("products", [])
            area_start = area.get("start")
            area_end = area.get("end")
            
            protoclusters = area.get("protoclusters", {})
            core_start = None
            core_end = None
            product_from_proto = None
            
            if protoclusters:
                first_proto = list(protoclusters.values())[0]
                core_start = first_proto.get("core_start")
                core_end = first_proto.get("core_end")
                product_from_proto = first_proto.get("product")
            
            # Region-specific gene function counts
            function_counts = Counter()
            if area_start is not None and area_end is not None:
                for gene_id, (g_start, g_end) in gene_locations.items():
                    # Check if gene overlaps with the BGC region
                    if g_start < area_end and g_end > area_start:
                        func = gene_to_function.get(gene_id)
                        if func:
                            function_counts[func] += 1
            
            results.append({
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
                "func_other": function_counts.get("other", 0)
            })

# Write CSV
fieldnames = [
    "genome_id", "record_id", "products", "start", "end",
    "core_start", "core_end", "knowncluster_hits",
    "func_transport", "func_regulatory", "func_biosynthetic",
    "func_biosynthetic_additional", "func_other"
]

with open(output_file, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"\nExtracted {len(results)} BGC regions from {len(json_files)} genomes")
print(f"Saved to: {output_file}")
