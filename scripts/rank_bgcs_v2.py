#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

IN_CSV = Path("results/04_bgc_metadata/bgc_with_metadata.csv")
OUT_DIR = Path("results/05_ranking")
OUT_DIR.mkdir(parents=True, exist_ok=True)

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
    return max(PRODUCT_WEIGHTS.get(p, 1.0) for p in parts)

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
    if gtype == "Isolate":
        score += 1.0
    elif gtype == "MAG":
        score += 0.5
    return score

def function_score(row):
    raw = (
        0.4 * row.get("func_transport", 0)
        + 0.3 * row.get("func_regulatory", 0)
        + 0.2 * row.get("func_biosynthetic", 0)
        + 0.2 * row.get("func_biosynthetic_additional", 0)
    )
    return min(float(raw), 3.0)

def has_class(products, name):
    if pd.isna(products):
        return False
    parts = [p.strip() for p in str(products).split(";")]
    return name in parts

df = pd.read_csv(IN_CSV)
df["score_novelty"] = df["knowncluster_hits"].apply(novelty_score)
df["score_product"] = df["products"].apply(product_score)
df["score_quality"] = df.apply(quality_score, axis=1)
df["score_function"] = df.apply(function_score, axis=1)
df["priority_score"] = (
    df["score_novelty"] + df["score_product"] + df["score_quality"] + df["score_function"]
)
df["rank_global"] = df["priority_score"].rank(ascending=False, method="min").astype(int)
df["rank_in_host"] = df.groupby("Host")["priority_score"].rank(ascending=False, method="min").astype(int)
df = df.sort_values(["Host", "priority_score"], ascending=[True, False])

df.to_csv(OUT_DIR / "bgc_ranked_v2.csv", index=False)

top25 = df[df["rank_in_host"] <= 25].copy()
top25.to_csv(OUT_DIR / "bgc_top25_per_host.csv", index=False)

for cls in ["hserlactone", "siderophore", "RiPP-like"]:
    sub = df[df["products"].apply(lambda x: has_class(x, cls))].copy()
    sub = sub.sort_values("priority_score", ascending=False)
    sub.head(50).to_csv(OUT_DIR / f"top50_{cls.replace('-', '_')}.csv", index=False)
    print(f"{cls}: {len(sub)} BGCs; top host counts:")
    print(sub.head(50)["Host"].value_counts())
    print()

print("Score min/median/max:", df["priority_score"].min(), df["priority_score"].median(), df["priority_score"].max())
print("\nTop 10 per host")
print(
    top25.groupby("Host")
    .head(10)[["Host", "genome_id", "Genome type", "products", "knowncluster_hits", "priority_score", "rank_in_host"]]
    .to_string(index=False)
)
