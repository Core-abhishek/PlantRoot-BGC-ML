#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

OUT = Path("results/05_ranking")
files = {
    "top25_per_host": OUT / "bgc_top25_per_host.csv",
    "hserlactone": OUT / "top50_hserlactone.csv",
    "siderophore": OUT / "top50_siderophore.csv",
    "RiPP-like": OUT / "top50_RiPP_like.csv",
}

cols = ["Host", "Phylum (gtdb)", "Family (gtdb)", "Genus (gtdb)", "Species (gtdb)"]

def show_counts(df, label):
    print(f"\n===== {label} ({len(df)} rows) =====")
    print("\nPhylum")
    print(df["Phylum (gtdb)"].value_counts(dropna=False).head(10).to_string())
    print("\nFamily")
    print(df["Family (gtdb)"].value_counts(dropna=False).head(10).to_string())
    print("\nGenus")
    print(df["Genus (gtdb)"].value_counts(dropna=False).head(15).to_string())
    if "Host" in df.columns:
        print("\nGenus x Host")
        tab = (
            df.groupby(["Host", "Genus (gtdb)"], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values("n", ascending=False)
        )
        print(tab.head(20).to_string(index=False))
        tab.to_csv(OUT / f"taxonomy_{label.replace('-', '_')}.csv", index=False)

for label, path in files.items():
    df = pd.read_csv(path)
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"{path} missing columns: {missing}")
        continue
    show_counts(df, label)
