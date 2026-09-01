# Plant root BGC Machine Learning Project

# PlantRoot-BGC-ML

Rule-based ranking of biosynthetic gene clusters (BGCs) from crop-root bacteria using the public CRBC dataset.

This is an active analysis pipeline, not a finished prediction tool.

## Aim

Build a reusable workflow that:

1. Extracts BGC features from antiSMASH JSON outputs
2. Joins them to CRBC genome metadata
3. Ranks clusters by novelty, product class, genome quality, and gene-function counts

The current ranking is an inventory shortlist. It does not prove plant-microbe communication.

## Current stage

- BGC feature table: 48,352 regions from 6,519 CRBC genomes
- Metadata join: host, genome type, quality, GTDB taxonomy
- 56.9% of BGCs have knowncluster hits = 0
- First rule-based ranking is available
- Next: correct ranking bias, then host-wise and class-wise shortlists

## Data

Source: Crop Root Bacterial Collection (CRBC) public release.

Large files are not stored in this repository:

- CRBC_BGC_part_1.tar.gz
- CRBC_BGC_part_2.tar.gz
- antiSMASH JSON files
- full merged or ranked tables if they are large

Keep those files locally under raw/ and processed/.

## Repository layout

scripts/
  extract_bgc_summary_batch.py
  extract_first_200_and_summarize.py
  extract_all_bgcs_and_summarize.py
  merge_bgc_metadata.py
  rank_bgcs.py
results/
  03_bgc_extraction/
  04_bgc_metadata/
  05_ranking/
progress.md

## Scoring used now

priority_score = novelty + product + quality + function

- Novelty: KnownClusterBlast hit count
- Product: higher weight for hserlactone, siderophore, NRPS/RiPP-like classes
- Quality: high-quality isolate scored above MAG
- Function: SMCOG counts inside the BGC region

Known limitation: large hybrid isolate BGCs are over-ranked. Rice novelty is under-used in the global top 100.

## Run

python3 scripts/merge_bgc_metadata.py
python3 scripts/rank_bgcs.py

JSON extraction from the CRBC tar archives:

python3 scripts/extract_all_bgcs_and_summarize.py

Requires Python 3, pandas, and openpyxl for the metadata Excel file.

## What this does not claim

- No experimental validation
- No metabolite structure
- No expression evidence
- No finished ML model
- No proof that a high-scoring BGC mediates plant-bacteria communication

## Status

Work in progress. See progress.md for the daily record.
