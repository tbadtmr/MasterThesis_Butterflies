#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH -J pisummary
#SBATCH -o logs/pi/pi_summary_%j.out
#SBATCH -e logs/pi/pi_summary_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="${REPO_ROOT}/07_model_comparison"

cd "$WORK"

module load PDC
module load R/4.5.2-cpeGNU-26.03

Rscript "$REPO_ROOT/00-scripts/10_summarize_validation_pi.R"
