#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J poolcheck
#SBATCH -o poolcheck_%j.out
#SBATCH -e poolcheck_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="${REPO_ROOT}/07_model_comparison"

cd "$WORK"

module load PDC
module load R/4.5.2-cpeGNU-26.03

echo "========================================"
echo "SAMPLING POOL CHECK"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "========================================"

Rscript --version
Rscript "$REPO_ROOT/00-scripts/05_check_sampling_pools.R"

echo
echo "Finished: $(date)"
