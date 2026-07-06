#!/bin/bash

# Paths
BIN_MOTIFS=/path/to/nanomotif/SAMPLE/bin-motifs.tsv
PILEUP=/path/to/modkit/SAMPLE.sorted.bed.gz
ASSEMBLY=/path/to/polished_assembly/SAMPLE.fasta
OUTDIR=/path/to/nanomotif/SAMPLE
MOTIFS_FILE=$OUTDIR/motifs.txt
OUTPUT=$OUTDIR/motifs-scored-read-methylation.tsv
THREADS=16

mkdir -p "$OUTDIR"

# motif list from bin-motifs.tsv (columns: bin, motif, mod_position, mod_type)
if [[ ! -s "$MOTIFS_FILE" ]]; then
    awk -F'\t' 'NR>1 {print $2 "_" $4 "_" $3}' "$BIN_MOTIFS" | sort -u > "$MOTIFS_FILE"
fi

MOTIFS=$(tr '\n' ' ' < "$MOTIFS_FILE")

epimetheus methylation-pattern contig \
    --pileup "$PILEUP" \
    --assembly "$ASSEMBLY" \
    --output "$OUTPUT" \
    --threads "$THREADS" \
    --batch-size 1000 \
    --output-type weighted-mean \
    --motifs $MOTIFS