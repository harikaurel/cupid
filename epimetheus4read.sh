#!/usr/bin/env bash

set -euo pipefail

mkdir -p epimetheus

# 1) Build the epimetheus motif list from the Nanomotif motifs table.

if [[ ! -s motifs.txt ]]; then
    awk -F'\t' '
        NR==1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "motif")                          m = i
                if ($i == "mod_type")                       t = i
                if ($i == "mod_position" || $i == "position") p = i
            }
            if (!(m && t && p)) {
                print "ERROR: motifs.tsv must have motif, mod_type, mod_position columns" > "/dev/stderr"
                exit 1
            }
            next
        }
        { print $m "_" $t "_" $p }
    ' motifs.tsv | sort -u > motifs.txt
    echo "[epimetheus4read] wrote $(wc -l < motifs.txt) motifs to motifs.txt"
fi

# 2) Pull the MM/ML (base-modification) tags out of the filtered BAM as FASTQ.
samtools fastq -T MM,ML filtered.bam > epimetheus/reads_MM_ML.fastq

# 3) Score per-read methylation for each motif.
epimetheus methylation-pattern \
    read-fastq --input epimetheus/reads_MM_ML.fastq \
    --output epimetheus/epimetheus_main.tsv \
    --motifs $(paste -sd' ' motifs.txt) \
    --threads 30