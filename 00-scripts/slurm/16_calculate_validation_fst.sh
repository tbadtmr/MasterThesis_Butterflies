#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 01:30:00
#SBATCH -J valfst
#SBATCH -o logs/fst/FST_%A_%a.out
#SBATCH -e logs/fst/FST_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="${REPO_ROOT}/07_model_comparison"

cd "$WORK"

module load PDC
module load R/4.5.2-cpeGNU-26.03

Rscript \
    "$REPO_ROOT/00-scripts/16_calculate_validation_fst.R" \
    "$SLURM_ARRAY_TASK_ID"
