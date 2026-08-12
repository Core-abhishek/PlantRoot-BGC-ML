import pandas as pd
import numpy as np

# -------------------------------------------------
# Load data
# -------------------------------------------------
df = pd.read_excel("CRBC_metadata_update.xlsx")

# Clean column names for easier use
df = df.rename(columns={
    "GenomeID_standard": "GenomeID",
    "Genome type": "Type",
    "Quality level": "Quality_level",
    "Phylum (gtdb)": "Phylum",
    "Family (gtdb)": "Family",
    "Genus (gtdb)": "Genus",
    "Class (gtdb)": "Class",
    "Order (gtdb)": "Order"
})

print("="*80)
print("CRBC METADATA – COMPREHENSIVE DISTRIBUTION ANALYSIS")
print("="*80)

# -------------------------------------------------
# 1. Total numbers
# -------------------------------------------------
print("\n1. TOTAL NUMBERS")
print("-"*50)
print(f"Total genomes: {len(df)}")
print("\nGenome type:")
print(df["Type"].value_counts())
print("\nQuality level:")
print(df["Quality_level"].value_counts())
print("\nHost:")
print(df["Host"].value_counts())

# -------------------------------------------------
# 2. Per host: number of genomes, types, quality levels
# -------------------------------------------------
print("\n\n2. PER HOST SUMMARY (Type + Quality)")
print("-"*50)
host_type_qual = pd.crosstab(
    index=df["Host"],
    columns=[df["Type"], df["Quality_level"]],
    margins=True
)
print(host_type_qual)

# -------------------------------------------------
# 3. Completeness distributions
# -------------------------------------------------
print("\n\n3. COMPLETENESS DISTRIBUTIONS")
print("-"*50)

def completeness_bins(series):
    return pd.cut(series,
                  bins=[-np.inf, 90, 95, np.inf],
                  labels=["<90%", "90-95%", "≥95%"])

df["Completeness_bin"] = completeness_bins(df["Completeness"])

print("\nOverall Completeness:")
print(df["Completeness_bin"].value_counts().sort_index())

print("\nCompleteness per Host:")
print(pd.crosstab(df["Host"], df["Completeness_bin"], margins=True))

print("\nCompleteness per Host × Type:")
print(pd.crosstab([df["Host"], df["Type"]], df["Completeness_bin"], margins=True))

print("\nCompleteness per Host × Type × Quality_level:")
print(pd.crosstab([df["Host"], df["Type"], df["Quality_level"]], 
                  df["Completeness_bin"], margins=True))

# -------------------------------------------------
# 4. Contamination distributions
# -------------------------------------------------
print("\n\n4. CONTAMINATION DISTRIBUTIONS")
print("-"*50)

def contamination_bins(series):
    return pd.cut(series,
                  bins=[-np.inf, 1, 5, np.inf],
                  labels=["<1%", "1-5%", ">5%"])

df["Contamination_bin"] = contamination_bins(df["Contamination"])

print("\nOverall Contamination:")
print(df["Contamination_bin"].value_counts().sort_index())

print("\nContamination per Host:")
print(pd.crosstab(df["Host"], df["Contamination_bin"], margins=True))

print("\nContamination per Host × Type:")
print(pd.crosstab([df["Host"], df["Type"]], df["Contamination_bin"], margins=True))

print("\nContamination per Host × Type × Quality_level:")
print(pd.crosstab([df["Host"], df["Type"], df["Quality_level"]], 
                  df["Contamination_bin"], margins=True))

# -------------------------------------------------
# 5. Quality score, GC, Coding density (same style)
# -------------------------------------------------
print("\n\n5. QUALITY SCORE / GC / CODING DENSITY (summary stats)")
print("-"*50)

for col in ["Quality score", "GC", "Coding density"]:
    print(f"\n{col} – Overall:")
    print(df[col].describe().round(3))
    print(f"\n{col} by Host × Type:")
    print(df.groupby(["Host", "Type"])[col].agg(["count", "mean", "median", "min", "max"]).round(3))

# -------------------------------------------------
# 6. Taxonomy (Phylum, Family, Genus)
# -------------------------------------------------
print("\n\n6. TAXONOMY")
print("-"*50)

print("\nTop 15 Phyla:")
print(df["Phylum"].value_counts().head(15))

print("\nTop 10 Families:")
print(df["Family"].value_counts().head(10))

print("\nTop 10 Genera:")
print(df["Genus"].value_counts().head(10))

print("\nPhylum by Host:")
print(pd.crosstab(df["Host"], df["Phylum"], margins=True))

# -------------------------------------------------
# 7. Save all important tables
# -------------------------------------------------
host_type_qual.to_csv("01_host_type_quality.csv")
pd.crosstab(df["Host"], df["Completeness_bin"], margins=True).to_csv("02_completeness_by_host.csv")
pd.crosstab([df["Host"], df["Type"]], df["Completeness_bin"], margins=True).to_csv("03_completeness_host_type.csv")
pd.crosstab(df["Host"], df["Contamination_bin"], margins=True).to_csv("04_contamination_by_host.csv")
pd.crosstab([df["Host"], df["Type"]], df["Contamination_bin"], margins=True).to_csv("05_contamination_host_type.csv")
pd.crosstab(df["Host"], df["Phylum"], margins=True).to_csv("06_phylum_by_host.csv")

print("\n\n" + "="*80)
print("All key tables have been saved as CSV files.")
print("="*80)
