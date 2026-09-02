2026-08-10: Successfully created GitHub repository and made first push.

2026-08-11: Started downloading CRBC files (BGC_part_1 & 2). gff, metadata and annotation files were downloaded from 'https://www.cropmicrobiome.com/' 
2026-08-12: Day 2 – Identified smaller CRBC files and started download of the files [metadata and gff files downloaded. Protein, gene and BGC_part_2 are currently being downloaded].
2026-08-12 (Day 2): 
- Downloaded metadata + GFF + started protein/genes downloads
- Explored full CRBC metadata
- Wrote and tested comprehensive metadata summary script
- Generated host, quality, completeness, contamination and taxonomy distribution tables

2026-08-14: Downloaded BGC parts 1 & 2, GFF, KEGG, proteins, metadata. Gene file downloading.
2026-08-18: Genes file downloaded. Basic inventory of all CRBC data files completed. Inspected genes archive structure.
2026-08-18: Confirmed structure of genes, proteins and GFF archives. All contain one compressed file per genome (~6700 files each).
2026-08-18: Inspected antiSMASH JSON structure from one example genome (CRBC_G0207).
Key fields identified: products, start/end, core_start/core_end, protoclusters (dict).
Ready to design bulk extraction of BGC summary table.

## 2026-08-25
- Generalised BGC extraction script to handle multiple JSON files
- Added knowncluster_hits (novelty indicator)
- Implemented region-specific gene function counts (SMCOGs)
- Tested successfully on 3 genomes
- Ready for bulk JSON extraction next

## ## 2026-08-27
- Optimised script for extracting 200 jsons in one batch, and extracting the csv summary from it.
- Optimised script to run on both the BGC tar.gz files, and run in loop for all the json files across the tar.gz file in batches of 200.
- The script was ran for the entire BGC file of CRBC, and the output summary file was generated. 

## 2026-08-30
- Merged BGC summary with CRBC metadata on genome_id = GenomeID_standard
- 48,352 BGC rows, 6,519 genomes, 0 rows missing host
- 27,513 BGCs (56.9%) have knowncluster_hits = 0
- Host BGC counts: Rice 18303, Wheat 14434, Maize 9502, Alfalfa 6113
- Zero-hit %: Rice 67.5, Alfalfa 55.7, Wheat 53.2, Maize 42.7
- MAG novelty 68.1% vs Isolate 51.3%
- Next: first BGC ranking score

## 2026-09-01 (ranking)
- First rule-based ranking on bgc_with_metadata.csv
- Score = novelty + product weight + quality + gene-function counts
- Score range 2.7 to 22.9
- Top 100 biased to maize isolate NRPS/PKS hybrids (Maize 62, Wheat 13, Alfalfa 13, Rice 12)
- This is an inventory shortlist, not communication ranking
- Next: cap function score, rank within host, separate lists for hserlactone, siderophore, RiPP-like

## 2026-09-02
- Corrected ranking: function score capped at 3, ranks within host
- Max score fell from 22.9 to 12.0
- Per-host top lists now dominated by novel isolate hserlactone BGCs
- Class shortlists: hserlactone 1885, siderophore 2283, RiPP-like 4692
- Taxonomy of top lists: hserlactone = Proteobacteria (Bradyrhizobium, Allorhizobium); siderophore more mixed (Solirubrobacter, Streptomyces, Bacillus, Acidovorax)
- Next: check unique genomes in the hserlactone top 5




























