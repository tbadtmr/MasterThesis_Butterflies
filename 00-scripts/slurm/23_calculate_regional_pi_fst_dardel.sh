#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -J popPiFst
#SBATCH -o logs/popPiFst_%j.out
#SBATCH -e logs/popPiFst_%j.err

set -euo pipefail

: "${WORK_DIR:?Set WORK_DIR before submitting}"
: "${DESIGN_FILE:?Set DESIGN_FILE before submitting}"
: "${REPLICATE_INVENTORY:?Set REPLICATE_INVENTORY before submitting}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COORD_FILE="${COORD_FILE:-$HOME/tabea_work/03_model_input/coord_chr9.txt}"
OUT_DIR="${OUT_DIR:-$WORK_DIR/results/regional_genetics}"

python3 "$SCRIPT_DIR/../23_calculate_regional_pi_fst.py" \
    --work-dir "$WORK_DIR" \
    --design "$DESIGN_FILE" \
    --replicate-inventory "$REPLICATE_INVENTORY" \
    --coord-file "$COORD_FILE" \
    --out-dir "$OUT_DIR"
