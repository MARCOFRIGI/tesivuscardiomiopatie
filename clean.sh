#!/bin/bash

set -euo pipefail

BASE_DIR="/mnt/d/pipeline_completa/Fase2_gnomAd"
GNOMAD_DIR="$BASE_DIR/esoma_gnomad"
REGION_DIR="$BASE_DIR/regions"
OUT_DIR="$BASE_DIR/gnomad_AF"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

CPU=$(nproc)
JOBS=16

echo "=========================================="
echo " Estrazione gnomAD v4.1 pulita"
echo " CPU disponibili: $CPU"
echo " Jobs in uso: $JOBS"
echo "=========================================="

process_chr() {

    region_file="$1"
    chr=$(basename "$region_file" | sed 's/regions_chr//' | sed 's/.txt//')

    vcf="$GNOMAD_DIR/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz"
    out="$OUT_DIR/chr${chr}.tsv"
    log="$LOG_DIR/chr${chr}.log"

    if [ ! -f "$vcf" ]; then
        echo "VCF non trovato per chr$chr" > "$log"
        return
    fi

    bcftools query \
        -R "$region_file" \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_grpmax\n' \
        "$vcf" > "$out" 2> "$log"
}

export -f process_chr
export GNOMAD_DIR REGION_DIR OUT_DIR LOG_DIR

parallel -j $JOBS --bar process_chr ::: "$REGION_DIR"/regions_chr*.txt

echo "Unione file..."

cat "$OUT_DIR"/chr*.tsv > gnomad_AF_all.tsv

echo "=========================================="
echo "File finale creato: gnomad_AF_all.tsv"
echo "=========================================="