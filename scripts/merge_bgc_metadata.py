#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BGC_CSV = Path("results/03_bgc_extraction/all_bgc_summary.csv")
META_XLSX = Path("raw/zenodo/CRBC_metadata_update.xlsx")
OUT_DIR = Path("results/04_bgc_metadata")
OUT_DIR.mkdir(parents=True, exist_ok=True)

bgc = pd.read_csv(BGC_CSV)
meta = pd.read_excel(META_XLSX)

meta_keep = meta[[
    "GenomeID_standard",
    "Host",
    "Genome type",
    "Quality level",
    "Completeness",
    "Contamination",
    "Quality score",
    "Phylum (gtdb)",
    "Class (gtdb)",
    "Order (gtdb)",
    "Family (gtdb)",
    "Genus (gtdb)",
    "Species (gtdb)",
    "Root compartment",
]].copy()
meta_keep = meta_keep.rename(columns={"GenomeID_standard": "genome_id"})

merged = bgc.merge(meta_keep, on="genome_id", how="left")
merged["is_novel"] = (merged["knowncluster_hits"].fillna(0) == 0).astype(int)

merged.to_csv(OUT_DIR / "bgc_with_metadata.csv", index=False)

print("BGC rows:", len(bgc))
print("Merged rows:", len(merged))
print("Merged genomes:", merged["genome_id"].nunique())
print("Rows missing host:", merged["Host"].isna().sum())
print("Zero-hit BGCs:", merged["is_novel"].sum())
print("Zero-hit percent:", round(100 * merged["is_novel"].mean(), 1))

novel = merged[merged["is_novel"] == 1].copy()

def count_table(df, cols, name):
    tab = (
        df.groupby(cols, dropna=False)
        .size()
        .reset_index(name="n_bgc")
        .sort_values("n_bgc", ascending=False)
    )
    tab.to_csv(OUT_DIR / name, index=False)
    return tab

count_table(merged, ["Host"], "bgc_by_host.csv")
count_table(novel, ["Host"], "novel_bgc_by_host.csv")
count_table(merged, ["Host", "Genome type"], "bgc_by_host_genometype.csv")
count_table(novel, ["Host", "Genome type"], "novel_bgc_by_host_genometype.csv")
count_table(merged, ["products"], "bgc_by_product.csv")
count_table(novel, ["products"], "novel_bgc_by_product.csv")
count_table(merged, ["Phylum (gtdb)"], "bgc_by_phylum.csv")
count_table(novel, ["Phylum (gtdb)"], "novel_bgc_by_phylum.csv")
count_table(merged, ["Host", "products"], "bgc_by_host_product.csv")
count_table(novel, ["Host", "products"], "novel_bgc_by_host_product.csv")

print("\nAll BGC by host")
print(merged.groupby("Host", dropna=False).size().sort_values(ascending=False))
print("\nZero-hit BGC by host")
print(novel.groupby("Host", dropna=False).size().sort_values(ascending=False))
print("\nAll BGC by genome type")
print(merged.groupby("Genome type", dropna=False).size().sort_values(ascending=False))
print("\nZero-hit BGC by genome type")
print(novel.groupby("Genome type", dropna=False).size().sort_values(ascending=False))
print("\nTop 15 product types overall")
print(merged["products"].value_counts().head(15))
print("\nTop 15 product types among zero-hit BGCs")
print(novel["products"].value_counts().head(15))
print(f"\nSaved merged table: {OUT_DIR / 'bgc_with_metadata.csv'}")
