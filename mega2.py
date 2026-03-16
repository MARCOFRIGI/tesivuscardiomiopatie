#!/usr/bin/env python3
"""
CARDIOMYOPATHY VARIANT DATASET — PIPELINE v6.2 GOLD STANDARD
Fix Header ClinVar Submission + Filtro Proteomico Strict + AF_bin categoriale + QC Tesi + Export Light
"""

import polars as pl
import pandas as pd
import numpy as np
import time, sys

T0 = time.time()
print("=" * 60)
print("  CARDIOMYOPATHY PIPELINE v6.2 (Gold Standard Edition)")
print("=" * 60)

# ─────────────────────────────────────────────────────────────
# CONFIGURAZIONE
# ─────────────────────────────────────────────────────────────
GENES = {
    "ACTC1","ALPK3","BAG3","CSRP3","DES","DSG2","DSP","FHOD3","FLNC",
    "JPH2","LDB3","LMNA","LMOD2","MYBPC3","MYH7","MYL2","MYL3","NEXN",
    "NRAP","PKP2","PLN","PPA2","PPP1R13L","PRDM16","PRKAG2","RBM20",
    "SCN5A","TBX20","TNNC1","TNNI3","TNNT2","TPM1","VCL"
}
print(f"  Target genes: {len(GENES)}")

def make_key(chrom_expr, pos_expr, ref_expr, alt_expr):
    return pl.concat_str([
        pl.lit("chr") + chrom_expr.str.replace_all(r"^chr", ""),
        pl.lit("_"), pos_expr.cast(pl.Utf8),
        pl.lit("_"), ref_expr.str.to_uppercase(),
        pl.lit("_"), alt_expr.str.to_uppercase()
    ])

def find_col(cols, *candidates):
    for c in candidates:
        if c in cols: return c
    return None

# ─────────────────────────────────────────────────────────────
# 0. GTF REGIONI ESONICHE
# ─────────────────────────────────────────────────────────────
print("\n[0/6] Regioni esoniche da GTF (opzionale)...")
try:
    gtf = pd.read_csv(
        "Homo_sapiens.GRCh38.115.gtf", sep="\t", comment="#", header=None,
        names=["chr","source","feature","start","end","score","strand","frame","attr"],
        dtype={"chr": str}
    )
    gtf = gtf[gtf["feature"] == "exon"].copy()
    gtf["gene_name"] = gtf["attr"].str.extract(r'gene_name "([^"]+)"')
    gtf = gtf[gtf["gene_name"].isin(GENES)]
    gtf["chr"] = np.where(gtf["chr"].str.startswith("chr"), gtf["chr"], "chr" + gtf["chr"])
    gtf[["chr","start","end"]].sort_values(["chr","start"]).drop_duplicates().to_csv(
        "GENE_REGIONS_33genes.tsv", sep="\t", header=False, index=False
    )
    print(f"  OK → GENE_REGIONS_33genes.tsv ({len(gtf):,} esoni)")
except FileNotFoundError:
    print("  GTF non trovato — saltato")

# ─────────────────────────────────────────────────────────────
# 1. VARIANT_SUMMARY + FILTRO PROTEOMICO
# ─────────────────────────────────────────────────────────────
print("\n[1/6] ClinVar variant_summary + Filtro Proteomico...")
vs = (
    pl.scan_csv("variant_summary.txt", separator="\t", infer_schema_length=0, ignore_errors=True)
    .filter(
        (pl.col("Assembly") == "GRCh38") & 
        (pl.col("Type") == "single nucleotide variant") & 
        (pl.col("PositionVCF") != "-1")
    )
    .with_columns([
        pl.col("Chromosome").str.replace_all(r"^chr","").alias("chrom_raw"),
        pl.col("PositionVCF").cast(pl.Int64, strict=False).alias("POS"),
        pl.col("ReferenceAlleleVCF").str.to_uppercase().alias("REF"),
        pl.col("AlternateAlleleVCF").str.to_uppercase().alias("ALT"),
        pl.col("VariationID").cast(pl.Int64, strict=False),
        pl.col("Name")
    ])
    .filter(
        pl.col("POS").is_not_null() & 
        (pl.col("REF").str.len_bytes() == 1) & 
        (pl.col("ALT").str.len_bytes() == 1) &
        pl.col("REF").str.contains(r"^[ACGT]$") &
        pl.col("ALT").str.contains(r"^[ACGT]$")
    )
    # FILTRO PROTEOMICO STRICT (Mantiene solo le varianti con impatto sull'amminoacido)
    .with_columns([
        pl.col("Name").str.extract(r'p\.[a-zA-Z]+(\d+)').cast(pl.Int64, strict=False).alias("aa_pos")
    ])
    .filter(pl.col("aa_pos").is_not_null()) # SCARTA TUTTE LE NON-MISSENSE/NON-CODING!
    .with_columns(pl.col("GeneSymbol").str.split(";").alias("_gl"))
    .explode("_gl")
    .with_columns(pl.col("_gl").str.strip_chars().alias("Gene"))
    .drop("_gl")
    .filter(pl.col("Gene").is_in(GENES))
    .with_columns(
        make_key(pl.col("chrom_raw"), pl.col("POS"), pl.col("REF"), pl.col("ALT")).alias("key")
    )
    .collect()
)

vs = vs.unique(subset=["key"], keep="first")
print(f"  Varianti SNV missense/protein-altering uniche: {len(vs):,}")
keys_lazy = vs.select("key").unique().lazy()

# ─────────────────────────────────────────────────────────────
# 2. FLAG CLINICI & INDICI DELLA TESI (PS, CI, PR)
# ─────────────────────────────────────────────────────────────
print("\n[2/6] Flag clinici ACMG e Indici (PS, CI, PR)...")

vs = vs.with_columns(
    pl.col("ClinicalSignificance").str.to_lowercase().alias("_cs")
).with_columns([
    (pl.col("_cs").str.contains("pathogenic") & ~pl.col("_cs").str.contains("likely pathogenic")).alias("is_pathogenic"),
    pl.col("_cs").str.contains("likely pathogenic").alias("is_likely_pathogenic"),
    pl.col("_cs").str.contains("uncertain significance").alias("is_vus"),
    (pl.col("_cs").str.contains("benign") & ~pl.col("_cs").str.contains("likely benign")).alias("is_benign"),
    pl.col("_cs").str.contains("likely benign").alias("is_likely_benign"),
    pl.col("_cs").str.contains("conflict").alias("is_conflicting"),
]).drop("_cs")

scv_col = find_col(vs.columns, "SCVsForAggregateGermlineClassification", "SCVsForAlleleGermlineClassification", "SCV")
if scv_col:
    vs = vs.with_columns(
        pl.when(pl.col(scv_col).is_null() | pl.col(scv_col).is_in(["-","."])).then(0)
          .otherwise(pl.col(scv_col).str.split("|").list.len()).alias("N_SCV")
    )
else:
    vs = vs.with_columns(pl.lit(0).alias("N_SCV"))

# FIX: Lettura ULTRA-STRICT dell'header di submission_summary.txt
_sub_header_idx = None
_sub_header_cols = None
try:
    with open("submission_summary.txt", "r", encoding="utf-8") as _fh:
        for _i, _line in enumerate(_fh):
            if _line.startswith("#VariationID") and "\t" in _line and "ClinicalSignificance" in _line:
                _sub_header_idx = _i
                _sub_header_cols = _line.lstrip("#").strip().split("\t")
                break
            if _i > 200:
                break
except FileNotFoundError:
    pass

if _sub_header_idx is None:
    print("  ⚠️ ERRORE: Header o file submission_summary.txt non trovato!")
    print("  Verranno impostati a NULL.")
    vs = vs.with_columns([
        pl.lit(None).cast(pl.Int64).alias("P_subs"), pl.lit(None).cast(pl.Int64).alias("B_subs"), 
        pl.lit(None).cast(pl.Int64).alias("Total_Committed"),
        pl.lit(None).cast(pl.Float64).alias("PR"), pl.lit(None).cast(pl.Float64).alias("PS"), pl.lit(None).cast(pl.Float64).alias("CI")
    ])
else:
    sub = pl.scan_csv(
        "submission_summary.txt", separator="\t", has_header=False,
        skip_rows=_sub_header_idx + 1, new_columns=_sub_header_cols, infer_schema_length=0, ignore_errors=True
    )
    
    sub_agg = (
        sub.with_columns(pl.col("VariationID").cast(pl.Int64, strict=False))
        .join(vs.select("VariationID").lazy(), on="VariationID", how="inner")
        .with_columns(pl.col("ClinicalSignificance").str.to_lowercase().alias("sig"))
        .with_columns([
            (pl.col("sig").str.contains("pathogenic") & ~pl.col("sig").str.contains("benign")).cast(pl.Int64).alias("is_P"),
            (pl.col("sig").str.contains("benign") & ~pl.col("sig").str.contains("pathogenic")).cast(pl.Int64).alias("is_B"),
        ])
        .group_by("VariationID")
        .agg([
            pl.col("is_P").sum().alias("P_subs"),
            pl.col("is_B").sum().alias("B_subs"),
        ])
        .with_columns((pl.col("P_subs") + pl.col("B_subs")).alias("Total_Committed"))
        .with_columns([
            pl.when(pl.col("Total_Committed") > 0).then(pl.col("P_subs") / pl.col("Total_Committed")).otherwise(None).alias("PR"),
            pl.when(pl.col("Total_Committed") > 0).then((pl.col("P_subs") - pl.col("B_subs")) / pl.col("Total_Committed")).otherwise(None).alias("PS")
        ])
        .with_columns([
            pl.when(pl.col("Total_Committed") == 0).then(None)
              .when((pl.col("PR") == 0.0) | (pl.col("PR") == 1.0)).then(0.0)
              .otherwise( - (pl.col("PR") * pl.col("PR").log(2)) - ((1.0 - pl.col("PR")) * (1.0 - pl.col("PR")).log(2)) )
              .alias("CI")
        ])
        .collect()
    )
    vs = vs.join(sub_agg, on="VariationID", how="left")
    print(f"  OK → Indici della tesi calcolati per {vs['PS'].drop_nulls().shape[0]:,} varianti.")

# ─────────────────────────────────────────────────────────────
# 3. GNOMAD & CATEGORIAL BINNING (AF_bin)
# ─────────────────────────────────────────────────────────────
print("\n[3/6] gnomAD AF_grpmax...")
try:
    gn_schema_cols = pl.scan_csv("UNIVERSE_gnomad_AF_full.tsv", separator="\t", infer_schema_length=10).collect_schema().names()
    chrom_gn  = find_col(gn_schema_cols, "CHROM","#CHROM","chrom")
    filter_gn = find_col(gn_schema_cols, "FILTER","filter")
    af_col    = find_col(gn_schema_cols, "AF","af")

    gn_expr = pl.scan_csv("UNIVERSE_gnomad_AF_full.tsv", separator="\t", infer_schema_length=0)
    if filter_gn:
        gn_expr = gn_expr.filter(pl.col(filter_gn) == "PASS")

    gn = (
        gn_expr
        .with_columns([
            pl.col(chrom_gn).str.replace_all(r"^chr","").alias("chrom_raw"),
            pl.col("POS").cast(pl.Int64, strict=False).alias("POS"),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
            pl.col("AF_grpmax").cast(pl.Float64, strict=False),
        ])
        .filter((pl.col("REF").str.len_bytes()==1) & (pl.col("ALT").str.len_bytes()==1))
        .with_columns(make_key(pl.col("chrom_raw"), pl.col("POS"), pl.col("REF"), pl.col("ALT")).alias("key"))
        .join(keys_lazy, on="key", how="inner")
        .select(["key","AF_grpmax"] + ([af_col] if af_col and af_col != "AF_grpmax" else []))
        .collect()
        .unique(subset=["key"])
    )

    vs = vs.join(gn, on="key", how="left")
except Exception as e:
    print(f"  ⚠️ ERRORE in gnomAD: {e}")
    vs = vs.with_columns(pl.lit(None).cast(pl.Float64).alias("AF_grpmax"))

vs = vs.with_columns([
    pl.when(pl.col("AF_grpmax").is_not_null())
      .then(pl.col("AF_grpmax").log(base=10.0))
      .otherwise(None)
      .alias("log10_AF_grpmax"),
    pl.when(pl.col("AF_grpmax").is_null() | (pl.col("AF_grpmax") < 1e-4)).then(pl.lit("ultra_rare"))
      .when((pl.col("AF_grpmax") >= 1e-4) & (pl.col("AF_grpmax") < 1e-3)).then(pl.lit("rare"))
      .otherwise(pl.lit("low_freq"))
      .alias("AF_bin")
])
print(f"  AF_grpmax matched: {vs['AF_grpmax'].drop_nulls().shape[0]:,} / {len(vs):,}")

# ─────────────────────────────────────────────────────────────
# 4. REVEL
# ─────────────────────────────────────────────────────────────
print("\n[4/6] REVEL scores...")
try:
    revel = (
        pl.scan_csv("revel_with_transcript_ids.tsv", separator=",",
                    infer_schema_length=0, null_values=".",
                    schema_overrides={"REVEL": pl.Float64, "grch38_pos": pl.Utf8})
        .with_columns([
            pl.col("chr").str.replace_all(r"^chr","").alias("chrom_raw"),
            pl.col("grch38_pos").cast(pl.Int64, strict=False).alias("POS"),
            pl.col("ref").str.to_uppercase().alias("REF"),
            pl.col("alt").str.to_uppercase().alias("ALT"),
        ])
        .with_columns(make_key(pl.col("chrom_raw"), pl.col("POS"), pl.col("REF"), pl.col("ALT")).alias("key"))
        .join(keys_lazy, on="key", how="inner")
        .group_by("key")
        .agg(pl.col("REVEL").max().alias("REVEL_max"))
        .collect()
    )

    vs = vs.join(revel, on="key", how="left")
except Exception as e:
    print(f"  ⚠️ ERRORE in REVEL: {e}")
    vs = vs.with_columns(pl.lit(None).cast(pl.Float64).alias("REVEL_max"))
    
print(f"  REVEL matched: {vs['REVEL_max'].drop_nulls().shape[0]:,} / {len(vs):,}")

# ─────────────────────────────────────────────────────────────
# 5. ALPHAMISSENSE (Con Soft-Fail per Robustezza)
# ─────────────────────────────────────────────────────────────
print("\n[5/6] AlphaMissense...")
_header_idx = None
_header_cols = None

try:
    with open("AlphaMissense_hg38.tsv", "r", encoding="utf-8") as _fh:
        for _i, _line in enumerate(_fh):
            _line_s = _line.rstrip("\n")
            if "\t" in _line_s and "am_pathogenicity" in _line_s:
                _header_idx = _i
                _header_cols = _line_s.lstrip("#").strip().split("\t")
                break
            if _i > 100:
                break
except FileNotFoundError:
    pass 

if _header_idx is None:
    print("  ⚠️ ERRORE: Header o file AlphaMissense non trovato. Verranno generati valori NULL.")
    alpha = pl.DataFrame({"key": [], "am_pathogenicity": [], "am_class": []}, 
                         schema={"key": pl.Utf8, "am_pathogenicity": pl.Float64, "am_class": pl.Utf8})
else:
    alpha = (
        pl.scan_csv(
            "AlphaMissense_hg38.tsv", separator="\t", has_header=False,
            skip_rows=_header_idx + 1, new_columns=_header_cols, infer_schema_length=0, null_values="."
        )
        .with_columns([
            pl.col("CHROM").str.replace_all(r"^chr", "").alias("chrom_raw"),
            pl.col("POS").cast(pl.Int64, strict=False).alias("POS"),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
            pl.col("am_pathogenicity").cast(pl.Float64, strict=False),
        ])
        .filter(pl.col("POS").is_not_null())
        .with_columns(make_key(pl.col("chrom_raw"), pl.col("POS"), pl.col("REF"), pl.col("ALT")).alias("key"))
        .join(keys_lazy, on="key", how="inner")
        .group_by("key")
        .agg([
            pl.col("am_pathogenicity").max(),
            pl.col("am_class").first()
        ])
        .collect()
    )

vs = vs.join(alpha, on="key", how="left")
print(f"  AlphaMissense matched: {vs['am_pathogenicity'].drop_nulls().shape[0]:,} / {len(vs):,}")

# ─────────────────────────────────────────────────────────────
# QC REPORT FINALE & EXPORT
# ─────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("  QC REPORT FINALE")
print("="*60)
print(f"  Varianti totali (solo Missense/Protein-altering): {len(vs):,}")
print(f"  PS / PR / CI disponibili:    {vs['PS'].drop_nulls().shape[0]:,}")
print(f"  AF_grpmax disponibile:       {vs['AF_grpmax'].drop_nulls().shape[0]:,}")
print(f"  REVEL disponibile:           {vs['REVEL_max'].drop_nulls().shape[0]:,}")
print(f"  AlphaMissense score:         {vs['am_pathogenicity'].drop_nulls().shape[0]:,}")

# Focus Tesi
vus_n = vs.filter(pl.col("is_vus")).height
print(f"\n  VUS focus tesi:              {vus_n:,} ({vus_n/len(vs)*100:.1f}%)")
if "PS" in vs.columns:
    ps_high_n = vs.filter(pl.col("PS") > 0.5).height
    print(f"  PS > 0.5 (high-risk):        {ps_high_n:,}")

print("\n[6/6] Export in corso...")

# 1. EXPORT DATASET COMPLETO
vs_pd = vs.to_pandas()
vs_pd.to_excel("cardiomyopathy_variants_FINAL_v6.xlsx", index=False)
vs.write_csv("cardiomyopathy_variants_FINAL_v6.csv")
print(f"  Dataset completo esportato -> cardiomyopathy_variants_FINAL_v6.xlsx/csv")

# 2. CREAZIONE ED EXPORT DATASET CONSULTABILE (LIGHT)
light_cols = [
    "Gene", "Name", "chrom_raw", "POS", "REF", "ALT", "RS# (dbSNP)", 
    "ClinicalSignificance", "N_SCV", "PS", "CI", 
    "AF_grpmax", "REVEL_max", "am_pathogenicity", "am_class"
]
# Manteniamo solo le colonne che esistono effettivamente nel dataframe per evitare errori
existing_light_cols = [c for c in light_cols if c in vs.columns]

# Selezioniamo le colonne e ordiniamo per Gene (A-Z) e REVEL (decrescente, in modo da avere le più gravi in alto)
vs_light = vs.select(existing_light_cols).sort(["Gene", "REVEL_max"], descending=[False, True], nulls_last=True)

vs_light_pd = vs_light.to_pandas()
vs_light_pd.to_excel("VUS_Consultabili_Score.xlsx", index=False)
vs_light.write_csv("VUS_Consultabili_Score.csv")
print(f"  Dataset consultabile esportato -> VUS_Consultabili_Score.xlsx/csv")

print(f"\n  Tempo totale: {round(time.time()-T0,1)} sec")