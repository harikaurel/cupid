#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Nanomotif: motif discovery + motif–contig scoring
# ============================================================

NANOMOTIF_BIN="/path/to/nanomotif"

ASSEMBLY_FASTA="/path/to/assembly.fasta"
PILEUP_BED="/path/to/pileup.bed"
CONTIG_BIN_TSV="/path/to/contig2bin.tsv"

OUT_DIR="/path/to/nanomotif_output"

THREADS=20
VALID_COVERAGE_THRESHOLD=1
MIN_MOTIF_SCORE=0.5
MIN_MOTIFS_BIN=1
THRESHOLD_VALID_COVERAGE=1

mkdir -p "$OUT_DIR"

# 1) Motif discovery
"$NANOMOTIF_BIN" motif_discovery \
  "$ASSEMBLY_FASTA" \
  "$PILEUP_BED" \
  -c "$CONTIG_BIN_TSV" \
  --out "$OUT_DIR" \
  -t "$THREADS" \
  --threshold_valid_coverage "$VALID_COVERAGE_THRESHOLD" \
    --min_motifs_bin "$MIN_MOTIFS_BIN" \
  --threshold_valid_coverage "$THRESHOLD_VALID_COVERAGE" \
  --min_motif_score "$MIN_MOTIF_SCORE" \
