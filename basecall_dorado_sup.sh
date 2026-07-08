#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Dorado basecalling (GPU)
# ============================================================

# ------------------------------------------------------------
# Input / output (edit as needed)
# ------------------------------------------------------------
READS_DIR="/path/to/pod5_dir"
OUTPUT_BAM="/path/to/output/trimmed.bam"
KIT_NAME="SQK-RBK114-24"

# ------------------------------------------------------------
# Tool paths (edit as needed)
# ------------------------------------------------------------
DORADO_BIN="/path/to/dorado/bin/dorado"
MODELS_DIR="/path/to/dorado/models"

# ------------------------------------------------------------
# Basecalling parameters
# ------------------------------------------------------------
MODEL_NAME="sup,6mA,4mC_5mC"
TRIM_MODE="all"
EMIT_MOVES=1
RECURSIVE=1

# ------------------------------------------------------------
# Prepare output directory
# ------------------------------------------------------------
mkdir -p "$(dirname "$OUTPUT_BAM")"

# ------------------------------------------------------------
# Build command
# ------------------------------------------------------------
CMD=("$DORADO_BIN" basecaller "$MODEL_NAME")

[[ "$RECURSIVE" == "1" ]] && CMD+=(-r)

CMD+=(
  "$READS_DIR"
  --models-directory "$MODELS_DIR"
  --kit-name "$KIT_NAME"
  --verbose
  --trim "$TRIM_MODE"
)

[[ "$EMIT_MOVES" == "1" ]] && CMD+=(--emit-moves)

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------
"${CMD[@]}" > "$OUTPUT_BAM"

"$DORADO_BIN" demux \
  --output-dir "$DEMUX_DIR" \
  --no-classify \
  "$OUTPUT_BAM"
