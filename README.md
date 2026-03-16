Pipeline di "Analisi delle VUS nelle Cardiomiopatie"
Tesi: "Oltre la Classificazione Standard: Un Approccio Integrato alla Valutazione delle Varianti di Significato Incerto nelle Cardiomiopatie" Autore: Marco Frigione Istituzione: Università Vita-Salute San Raffaele, Milano Supervisori: Prof. Giuseppe Banfi, Prof.ssa Chiara Di Resta, Prof. Giovanni Peretto

Panoramica
Questa repository contiene la pipeline di analisi utilizzata nella mia tesi per studiare 42.726 varianti missenso in 33 geni associati alle cardiomiopatie. La pipeline esegue un clustering GMM non supervisionato nello spazio dei punteggi AlphaMissense x REVEL e valida il gradiente di rischio risultante attraverso evidenze biologiche, popolazionali e temporali indipendenti.

La pipeline è interamente in locale. Durante l'esecuzione non vengono usate API. A parità di file di input, tutti gli output sono riproducibili.

File in questa repository
scripts/mega2.py — costruzione del dataset (passaggi 1-5) scripts/annotazione_unificata.py — annotazione strutturale (passaggi 6-8) scripts/statistica_finale.py — clustering e analisi statistica (16 sezioni) scripts/validazione_temporale_clinvar3.py — validazione temporale su ClinVar 2015-2025 scripts/identikit.py — validazione dell'identikit gene-specifico (fingerprint) bash/extract_gnomad_parallel.sh — estrazione delle frequenze alleliche da gnomAD bash/clean.sh — estrazione da gnomAD (versione pulita) data/README_data.md — fonti dei dati e istruzioni per il download requirements.txt — dipendenze Python

Pipeline — 8 passaggi
Passaggio 0 (bash): Estrae le coordinate esoniche per i 33 geni target dal file GTF di Ensembl (GRCh38 release 115) utilizzando extract_gnomad_parallel.sh. Utilizza bcftools con il flag -R per estrarre i dati di AF (frequenza allelica) dai file VCF di gnomAD v4.1, parallelizzati con GNU parallel. Output: UNIVERSE_gnomad_AF_full.tsv.

Passaggio 1 (mega2.py): Carica variant_summary.txt da ClinVar (Febbraio 2025). Applica filtri sequenziali: assemblaggio GRCh38, tipo SNV, posizione VCF valida, singolo nucleotide REF/ALT, filtro protein-altering tramite regex sul campo HGVS Name, filtro sui geni (33 target), deduplicazione sulla chiave genomica chr_POS_REF_ALT.

Passaggio 2 (mega2.py): Deriva 6 flag booleani da ClinicalSignificance con logica mutuamente esclusiva. Le interpretazioni discordanti (conflicting) sovrascrivono tutti gli altri flag. Calcola gli indici di sottomissione PS, PR, CI da submission_summary.txt.

Passaggio 3 (mega2.py): Left-join di UNIVERSE_gnomad_AF_full.tsv sulla chiave genomica. Mantiene AF_grpmax (frequenza allelica massima tra le popolazioni gnomAD). Calcola log10_AF_grpmax per le analisi a valle.

Passaggio 4 (mega2.py): Join di revel_with_transcript_ids.tsv (REVEL v1.3) sulla chiave genomica. Mantiene il punteggio REVEL massimo tra i trascritti (REVEL_max).

Passaggio 5 (mega2.py): Parsing di AlphaMissense_hg38.tsv (Ottobre 2023) con rilevamento robusto dell'intestazione. Mantiene la am_pathogenicity massima tra i trascritti. Output: cardiomyopathy_variants_FINAL_v6.xlsx.

Passaggi 6-7 (annotazione_unificata.py): Carica il TSV massivo UniProt Swiss-Prot (Aprile 2025, solo proteine revisionate). Mappa il gene all'ID UniProt per i 33 target. Estrae gli intervalli delle feature posizionali (domain, coiled_coil, repeat, region, motif, transmembrane, binding_site, active_site, ptm) tramite parsing regex dei campi FT di UniProt. Assegna flag binari per variante in base alla sovrapposizione della posizione amminoacidica.

Passaggio 8 (annotazione_unificata.py): Scansiona protein2ipr.dat.gz (InterPro release 100.0) filtrando per gli ID UniProt target. Deduplica le voci per (ipr_id, start, end). Assegna in_interpro, is_structural_domain, is_regulatory_region, is_family tramite corrispondenza posizionale. Output: cardiomyopathy_variants_with_interpro_local.xlsx.

Analisi statistica (statistica_finale.py)
Lo script principale esegue 16 sezioni sequenziali sul dataset annotato:

Caricamento del dataset e fix ClinVar (conflicting = priorità massima)
Clustering GMM: selezione del modello BIC k=1-8 con criterio del gomito (elbow) al 10%, normalizzazione StandardScaler
Riordinamento dei cluster per mediana AlphaMissense crescente (C0=basso rischio, C3=alto rischio)
Label propagation per 14.028 varianti prive di entrambi i punteggi (centroide più vicino su aa_pos_norm)
Gradiente di patogenicità: correlazione di Spearman, OR di Fisher (HR vs C0)
Validazione della frequenza allelica: Kruskal-Wallis, Mann-Whitney, Spearman cluster x log10AF
Distribuzione delle VUS: OR di Fisher per la sovrarappresentazione del cluster intermedio
Arricchimento dei domini proteici: test esatto di Fisher + Benjamini-Hochberg FDR
Arricchimento delle classi funzionali (sarcomeriche, desmosomiali, citoscheletriche, canali ionici)
AUC-ROC gene-specifica (punteggio composito sqrt(AM x REVEL) vs AM vs REVEL)
Analisi della banda IEV (AF da 1e-4 a 1e-2)
Robustezza: bootstrap n=100, stabilità del seed n=10, cross-validation n=20 (ARI)
Identikit patogenetico gene-specifico (fingerprint): intervalli percentili 10-90, domini arricchiti, soglia di AF
Scoring delle VUS candidate: corrispondenza dei criteri >= 80%, confronto con la soglia globale
Paesaggio di probabilità (Probability landscape): griglia 20x20 P(patogenetica) nello spazio AM x REVEL
Figure 1-8
Tutti i risultati numerici riportati nella tesi sono tracciabili nei file CSV presenti nella directory PIPELINE_RESULTS/. Il file master_metrics.csv contiene tutti i valori chiave in una singola riga verificabile.

Validazione temporale (validazione_temporale_clinvar3.py)
Validazione quasi-prospettica su 11 release annuali di ClinVar (2015-2025). Abbinamento tramite VariantID (2015-2016) e AlleleID (2017-2025, a seguito del cambio di formato di ClinVar). Per ogni anno: identifica le VUS di quella release riclassificate come P/LP nel dataset attuale, testa se l'appartenenza al cluster ad alto rischio predice la riclassificazione (OR di Fisher), valida gli 8 pattern biologici della tesi.

Input richiesto: variant_summary_YYYY.txt.gz per ogni anno dal 2015 al 2025.

Validazione dell'identikit gene-specifico (identikit.py)
Confronta tre metodi di prioritizzazione su tutte le 11 release annuali di ClinVar: l'identikit gene-specifico (corrispondenza dei criteri >= 80%), la soglia globale (AM >= 0.564 e REVEL >= 0.75 e AF < 1e-4) e l'appartenenza al cluster ad alto rischio (C3). Metriche per anno: sensibilità, valore predittivo positivo (PPV), OR di Fisher.

Fonti dei dati
Tutti i file di input devono essere scaricati in autonomia. Vedi data/README_data.md per le versioni esatte e gli URL.

variant_summary.txt e submission_summary.txt: FTP ClinVar (NCBI), Febbraio 2025 gnomad.exomes.v4.1.sites.chr*.vcf.bgz: browser gnomAD, versione 4.1 (GRCh38) Homo_sapiens.GRCh38.115.gtf: FTP Ensembl, release 115 revel_with_transcript_ids.tsv: sito web REVEL, versione 1.3 AlphaMissense_hg38.tsv: Google DeepMind Cloud Storage, Ottobre 2023 uniprotkb_*.tsv.gz: download massivo UniProt, Aprile 2025 (Swiss-Prot revisionate) protein2ipr.dat.gz e entry.list: FTP InterPro (EBI), release 100.0 variant_summary_YYYY.txt.gz: archivio FTP ClinVar, release annuali 2015-2025

Requisiti
Python >= 3.12 pandas >= 2.0 numpy >= 1.24 scipy >= 1.11 scikit-learn >= 1.3 matplotlib >= 3.7 seaborn >= 0.12 polars >= 0.18 openpyxl >= 3.1

Installa con: pip install -r requirements.txt

Utilità Bash: bcftools >= 1.17, GNU parallel
Fonti dei Dati

Tutti i file di input devono essere scaricati in autonomia prima di eseguire la pipeline.
Nessun file di dati è incluso in questa repository a causa delle dimensioni e dei vincoli di licenza.

File Necessari

ClinVar
File: variant_summary.txt, submission_summary.txt
Versione utilizzata: Febbraio 2025
Download:

gnomAD Exomes v4.1
File: gnomad.exomes.v4.1.sites.chr{1-22,X,Y}.vcf.bgz (+ file indice .tbi)
Versione utilizzata: gnomAD v4.1 (GRCh38)
Download:
Nota: Ogni file cromosomico pesa 5–30 GB. Vengono utilizzate solo le regioni esoniche dei 33 geni target — gli script di estrazione bash (Passaggio 0) usano bcftools con -R per estrarre solo le posizioni rilevanti.

Ensembl GTF (per estrazione regioni esoniche)
File: Homo_sapiens.GRCh38.115.gtf
Versione utilizzata: Ensembl release 115
Download:

REVEL
File: revel_with_transcript_ids.tsv
Versione utilizzata: v1.3
Download:

AlphaMissense
File: AlphaMissense_hg38.tsv
Versione utilizzata: Ottobre 2023
Download:
Nota: Richiede un account Google. Il file contiene i punteggi per tutte le possibili varianti missenso umane su hg38.

UniProt Swiss-Prot
File: uniprotkb_reviewed_*.tsv.gz (download massivo)
Versione utilizzata: Aprile 2025 (solo Swiss-Prot revisionate)
Download:

Vai su

Filtro: Reviewed (Swiss-Prot), Organism: Homo sapiens

Seleziona le colonne: Entry, Entry Name, Gene Names, Length, Involvement in disease, Domain [FT], Coiled coil, Repeat, Region, Motif, Transmembrane, Active site, Binding site, Modified residue, Protein families, Gene Ontology (molecular function), Gene Ontology (biological process), Subunit structure

Scarica come TSV (compresso)

InterPro
File: protein2ipr.dat.gz, entry.list
Versione utilizzata: Release 100.0
Download:

Nota: protein2ipr.dat.gz è ~20 GB non compresso. Lo script lo scansiona filtrando per i 33 ID UniProt target.

Archivio Storico ClinVar (per validazione temporale)
File: variant_summary_YYYY*.txt.gz per gli anni 2015–2025
Download:
Nota: Snapshot annuali disponibili dal 2015. La nomenclatura dei file varia in base all'anno.

Posizionamento dei File

Posiziona tutti i file nella directory di lavoro (working directory) insieme agli script, oppure aggiorna i percorsi all'inizio di ogni script. Gli script rilevano automaticamente i nomi dei file che corrispondono ai pattern attesi.

Requisiti di Spazio su Disco

Componente: Dimensione approssimativa

gnomAD v4.1 (tutti i cromosomi): ~400 GB

ClinVar variant_summary: ~200 MB

AlphaMissense hg38: ~6 GB

REVEL: ~6 GB

protein2ipr.dat.gz: ~6 GB compresso

UniProt TSV: ~500 MB
Totale: ~420 GB

Il dataset annotato finale (cardiomyopathy_variants_with_interpro_local.xlsx) è di ~50 MB.

Come citare
Frigione M. "Beyond Standard Classification: An Integrated Approach to the Evaluation of Variants of Uncertain
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
Data Sources

All input files must be downloaded independently before running the pipeline.
No data files are included in this repository due to size and licensing constraints.

Required Files

ClinVar
Files: variant_summary.txt, submission_summary.txt
Version used: February 2025
Download:


gnomAD Exomes v4.1
Files: gnomad.exomes.v4.1.sites.chr{1-22,X,Y}.vcf.bgz (+ .tbi index files)
Version used: gnomAD v4.1 (GRCh38)
Download: 
Note: Each chromosome file is 5–30 GB. Only exonic regions of the 33 target genes are used — the bash extraction scripts (Step 0) use bcftools with -R to extract only relevant positions.

Ensembl GTF (for exonic region extraction)
File: Homo_sapiens.GRCh38.115.gtf
Version used: Ensembl release 115
Download:

REVEL
File: revel_with_transcript_ids.tsv
Version used: v1.3
Download: 

AlphaMissense
File: AlphaMissense_hg38.tsv
Version used: October 2023
Download: 
Note: Requires Google account. File contains scores for all possible human missense variants on hg38.

UniProt Swiss-Prot
File: uniprotkb_reviewed_*.tsv.gz (bulk download)
Version used: April 2025 (Swiss-Prot reviewed only)
Download:

Go to 

Filter: Reviewed (Swiss-Prot), Organism: Homo sapiens

Select columns: Entry, Entry Name, Gene Names, Length, Involvement in disease, Domain [FT], Coiled coil, Repeat, Region, Motif, Transmembrane, Active site, Binding site, Modified residue, Protein families, Gene Ontology (molecular function), Gene Ontology (biological process), Subunit structure

Download as TSV (compressed)

InterPro
Files: protein2ipr.dat.gz, entry.list
Version used: Release 100.0
Download:


Note: protein2ipr.dat.gz is ~20 GB uncompressed. The script scans it filtering for the 33 target UniProt IDs.

ClinVar Historical Releases (for temporal validation)
Files: variant_summary_YYYY*.txt.gz for years 2015–2025
Download: 
Note: Annual snapshots available from 2015. File naming varies by year.

File Placement

Place all files in the working directory alongside the scripts, or update the paths at the top of each script. The scripts auto-detect filenames matching expected patterns.

Disk Space Requirements

Component: Approximate size

gnomAD v4.1 (all chromosomes): ~400 GB

ClinVar variant_summary: ~200 MB

AlphaMissense hg38: ~6 GB

REVEL: ~6 GB

protein2ipr.dat.gz: ~6 GB compressed

UniProt TSV: ~500 MB
Total: ~420 GB

The final annotated dataset (cardiomyopathy_variants_with_interpro_local.xlsx) is ~50 MB.

How to cite
Frigione M. "Beyond Standard Classification: An Integrated Approach to the Evaluation of Variants of Uncertain Significance in Cardiomyopathies." Master's Thesis, Università Vita-Salute San Raffaele, 2026.
