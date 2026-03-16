Cardiomyopathy VUS Analysis Pipeline
Thesis: "Beyond Standard Classification: An Integrated Approach to the Evaluation of Variants of Uncertain Significance in Cardiomyopathies" Author: Marco Frigione Institution: Università Vita-Salute San Raffaele, Milan Supervisors: Prof. Giuseppe Banfi, Prof.ssa Chiara Di Resta, Prof. Giovanni Peretto

Overview
This repository contains the analysis pipeline used in my thesis to study 42,726 missense variants across 33 cardiomyopathy-associated genes. The pipeline performs unsupervised GMM clustering in the AlphaMissense x REVEL score space and validates the resulting risk gradient through independent biological, populational, and temporal evidence.

The pipeline is fully local. No network calls are made during execution. Given the same input files, all outputs are reproducible (fixed random seeds throughout).

Files in this repository
scripts/mega2.py — dataset construction (steps 1-5) scripts/annotazione_unificata.py — structural annotation (steps 6-8) scripts/statistica_finale.py — clustering and statistical analysis (16 sections) scripts/validazione_temporale_clinvar3.py — temporal validation across ClinVar 2015-2025 scripts/identikit.py — gene-specific fingerprint validation bash/extract_gnomad_parallel.sh — gnomAD allele frequency extraction bash/clean.sh — gnomAD extraction (clean version) data/README_data.md — data sources and download instructions requirements.txt — Python dependencies

Pipeline — 8 steps
Step 0 (bash): Extract exonic coordinates for the 33 target genes from the Ensembl GTF file (GRCh38 release 115) using extract_gnomad_parallel.sh. Uses bcftools with -R flag to extract AF data from gnomAD v4.1 VCF files, parallelised with GNU parallel. Output: UNIVERSE_gnomad_AF_full.tsv.

Step 1 (mega2.py): Load variant_summary.txt from ClinVar (February 2025). Apply sequential filters: GRCh38 assembly, SNV type, valid VCF position, single nucleotide REF/ALT, protein-altering filter via regex on HGVS Name field, gene filter (33 targets), deduplication on genomic key chr_POS_REF_ALT.

Step 2 (mega2.py): Derive 6 boolean flags from ClinicalSignificance with mutually exclusive logic. Conflicting interpretations override all other flags. Calculate submission indices PS, PR, CI from submission_summary.txt.

Step 3 (mega2.py): Left-join UNIVERSE_gnomad_AF_full.tsv on genomic key. Retain AF_grpmax (maximum allele frequency across gnomAD populations). Compute log10_AF_grpmax for downstream analyses.

Step 4 (mega2.py): Join revel_with_transcript_ids.tsv (REVEL v1.3) on genomic key. Retain maximum REVEL score across transcripts (REVEL_max).

Step 5 (mega2.py): Parse AlphaMissense_hg38.tsv (October 2023) with robust header detection. Retain maximum am_pathogenicity across transcripts. Output: cardiomyopathy_variants_FINAL_v6.xlsx.

Steps 6-7 (annotazione_unificata.py): Load UniProt Swiss-Prot bulk TSV (April 2025, reviewed proteins only). Map gene to UniProt ID for 33 targets. Extract positional feature intervals (domain, coiled_coil, repeat, region, motif, transmembrane, binding_site, active_site, ptm) via regex parsing of UniProt FT fields. Assign binary flags per variant based on amino acid position overlap.

Step 8 (annotazione_unificata.py): Scan protein2ipr.dat.gz (InterPro release 100.0) filtering for target UniProt IDs. Deduplicate entries by (ipr_id, start, end). Assign in_interpro, is_structural_domain, is_regulatory_region, is_family via positional matching. Output: cardiomyopathy_variants_with_interpro_local.xlsx.

Statistical analysis (statistica_finale.py)
The main script runs 16 sequential sections on the annotated dataset:

Dataset loading and ClinVar fix (conflicting = maximum priority)
GMM clustering: BIC model selection k=1-8 with 10% elbow criterion, StandardScaler normalisation
Cluster reordering by ascending AlphaMissense median (C0=low risk, C3=high risk)
Label propagation for 14,028 variants missing both scores (nearest centroid on aa_pos_norm)
Pathogenicity gradient: Spearman correlation, Fisher OR (HR vs C0)
Allele frequency validation: Kruskal-Wallis, Mann-Whitney, Spearman cluster x log10AF
VUS distribution: Fisher OR for intermediate cluster overrepresentation
Protein domain enrichment: Fisher exact + Benjamini-Hochberg FDR
Functional class enrichment (sarcomeric, desmosomal, cytoskeletal, ion channel)
Gene-specific AUC-ROC (composite score sqrt(AM x REVEL) vs AM vs REVEL)
IEV band analysis (AF 1e-4 to 1e-2)
Robustness: bootstrap n=100, seed stability n=10, cross-validation n=20 (ARI)
Gene-specific pathogenic fingerprint: percentile 10-90 intervals, enriched domains, AF threshold
VUS candidate scoring: >= 80% criteria match, comparison with global threshold
Probability landscape: 20x20 grid P(pathogenic) in AM x REVEL space
Figures 1-8
All numerical results reported in the thesis trace back to CSV files in the PIPELINE_RESULTS/ directory. The file master_metrics.csv contains all key values as a single verifiable row.

Temporal validation (validazione_temporale_clinvar3.py)
Quasi-prospective validation across 11 annual ClinVar releases (2015-2025). Matching via VariantID (2015-2016) and AlleleID (2017-2025, following ClinVar format change). For each year: identifies VUS in that release reclassified as P/LP in the current dataset, tests whether high-risk cluster membership predicts reclassification (Fisher OR), validates 8 biological patterns from the thesis.

Required input: variant_summary_YYYY.txt.gz for each year from 2015 to 2025.

Gene-specific fingerprint validation (identikit.py)
Compares three prioritisation methods across all 11 annual ClinVar releases: the gene-specific identikit (>= 80% criteria match), the global threshold (AM >= 0.564 and REVEL >= 0.75 and AF < 1e-4), and high-risk cluster membership (C3). Metrics per year: sensitivity, positive predictive value, Fisher OR.

Data sources
All input files must be downloaded independently. See data/README_data.md for exact versions and URLs.

variant_summary.txt and submission_summary.txt: ClinVar FTP (NCBI), February 2025 gnomad.exomes.v4.1.sites.chr*.vcf.bgz: gnomAD browser, version 4.1 (GRCh38) Homo_sapiens.GRCh38.115.gtf: Ensembl FTP, release 115 revel_with_transcript_ids.tsv: REVEL website, version 1.3 AlphaMissense_hg38.tsv: Google DeepMind Cloud Storage, October 2023 uniprotkb_*.tsv.gz: UniProt bulk download, April 2025 (Swiss-Prot reviewed) protein2ipr.dat.gz and entry.list: InterPro FTP (EBI), release 100.0 variant_summary_YYYY.txt.gz: ClinVar FTP archive, annual releases 2015-2025

Requirements
Python >= 3.12 pandas >= 2.0 numpy >= 1.24 scipy >= 1.11 scikit-learn >= 1.3 matplotlib >= 3.7 seaborn >= 0.12 polars >= 0.18 openpyxl >= 3.1

Install with: pip install -r requirements.txt

Bash utilities: bcftools >= 1.17, GNU parallel

How to cite
Frigione M. "Beyond Standard Classification: An Integrated Approach to the Evaluation of Variants of Uncertain Significance in Cardiomyopathies." Master's Thesis, Università Vita-Salute San Raffaele, 2026.
