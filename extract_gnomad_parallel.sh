#!/bin/bash

# ==========================================
# ESTRAZIONE AF DA GNOMAD EXOMES v4.1
# Pipeline stabile con GNU parallel
# ==========================================

BASE_DIR="/mnt/d/pipeline_completa/Fase2_gnomAd"
GNOMAD_DIR="$BASE_DIR/esoma_gnomad"
REGION_DIR="$BASE_DIR/regions"
OUT_DIR="$BASE_DIR/gnomad_AF"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

echo "=========================================="
echo "Estrazione gnomAD v4.1 avviata"
echo "CPU disponibili: $(nproc)"
echo "Jobs in uso: 16"
echo "=========================================="

process_chr() {

    chr_file="$1"
    chr=$(basename "$chr_file" | sed 's/regions_chr//' | sed 's/.txt//')

    vcf_file="$GNOMAD_DIR/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz"
    output_file="$OUT_DIR/AF_chr${chr}.tsv"
    log_file="$LOG_DIR/log_chr${chr}.txt"

    echo "Starting chr${chr}..."

    if [ -f "$vcf_file" ]; then

        bcftools query \
            -R "$chr_file" \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_grpmax\n' \
            "$vcf_file" \
            > "$output_file" 2> "$log_file"

        echo "Completed chr${chr}"

    else
        echo "VCF non trovato per chr${chr}" > "$log_file"
    fi
}

export -f process_chr
export GNOMAD_DIR REGION_DIR OUT_DIR LOG_DIR

parallel --jobs 16 --bar process_chr ::: "$REGION_DIR"/regions_chr*.txt

echo "=========================================="
echo "Unione file..."
echo "=========================================="

cat "$OUT_DIR"/AF_chr*.tsv > "$BASE_DIR/gnomad_AF_all.tsv"

echo "File finale creato: gnomad_AF_all.tsv"
echo "Output directory per cromosoma: $OUT_DIR"
echo "=========================================="
echo "Estrazione completata"
echo "=========================================="