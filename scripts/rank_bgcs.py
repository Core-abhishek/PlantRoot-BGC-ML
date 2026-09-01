#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

IN_CSV = Path("results/04_bgc_metadata/bgc_with_metadata.csv")
OUT_DIR = Path("results/05_ranking")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Higher weight = more relevant to interaction / signalling / antagonism
PRODUCT_WEIGHTS = {
    "hserlactone": 3.0,
    "siderophore": 2.5,
    "NRPS": 2.0,
    "NRPS-like": 2.0,
    "RiPP-like": 2.0,
    "lassopeptide": 2.0,
    "thioamitides": 1.8,
    "arylpolyene": 1.5,
    "betalactone": 1.5,
    "T1PKS": 1.5,
    "T3PKS": 1.5,
    "terpene": 1.2,
    "RRE-containing": 1.2,
    "redox-cofactor": 1.0,
    "NAPAA": 1.0,
}

def product_score(products):
    if pd.isna(products) or str(products).strip() == "":
        return 1.0
    parts = [p.strip() for p in str(products).split(";")]
    scores = [PRODUCT_WEIGHTS.get(p, 1.0) for p in parts]
    return max(scores)

def novelty_score(hits):
    hits = 0 if pd.isna(hits) else int(hits)
    if hits == 0:
        return 3.0
    if hits <= 2:
        return 2.0
    if hits <= 5:
        return 1.0
    return 0.0

def quality_score(row):
    level = str(row.get("Quality level", "")).lower()
    gtype = str(row.get("Genome type", ""))
    score = 0.0
    if "high" in level:
        score += 2.0
    elif "medium" in level:
        score += 1.0
    else:
        score += 0.0
    if gtype == "Isolate":
        score += 1.0
    elif gtype == "MAG":
        score += 0.5
    return score

def function_score(row):
    return (
        0.4 * row.get("func_transport", 0)
        + 0.3 * row.get("func_regulatory", 0)
        + 0.2 * row.get("func_biosynthetic", 0)
        + 0.2 * row.get("func_biosynthetic_additional", 0)
    )

df = pd.read_csv(IN_CSV)
df["score_novelty"] = df["knowncluster_hits"].apply(novelty_score)
df["score_product"] = df["products"].apply(product_score)
df["score_quality"] = df.apply(quality_score, axis=1)
df["score_function"] = df.apply(function_score, axis=1)
df["priority_score"] = (
    df["score_novelty"]
    + df["score_product"]
    + df["score_quality"]
    + df["score_function"]
)
df = df.sort_values("priority_score", ascending=False)

df.to_csv(OUT_DIR / "bgc_ranked.csv", index=False)
df.head(100).to_csv(OUT_DIR / "bgc_top100.csv", index=False)

print("Ranked BGCs:", len(df))
print("Score min/median/max:", df["priority_score"].min(), df["priority_score"].median(), df["priority_score"].max())
print("\nTop 15")
print(
    df[
        ["genome_id", "Host", "Genome type", "products", "knowncluster_hits", "priority_score"]
    ].head(15).to_string(index=False)
)
print("\nTop 100 by host")
print(df.head(100)["Host"].value_counts())
print("\nTop 100 by product")
print(df.head(100)["products"].value_counts().head(15))
print(f"\nSaved: {OUT_DIR / 'bgc_ranked.csv'}")
print(f"Saved: {OUT_DIR / 'bgc_top100.csv'}")
