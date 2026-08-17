#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J valfroh_r
#SBATCH -o logs/froh/FROH_%A_%a_rate.out
#SBATCH -e logs/froh/FROH_%A_%a_rate.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="${REPO_ROOT}/07_model_comparison"

cd "$WORK"

ROH_RATE="${ROH_RATE:-1.3e-8}"
ROH_TAG="${ROH_TAG:-rate1p3e8}"
OUTROOT="results/froh_${ROH_TAG}"

TASK="$SLURM_ARRAY_TASK_ID"
TASKPAD=$(printf "%03d" "$TASK")

mkdir -p "${OUTROOT}/raw"
mkdir -p "${OUTROOT}/per_snapshot"
mkdir -p logs/froh

# ============================================================
# 1. GET TASK METADATA
# ============================================================

read LANDSCAPE MODEL REP YEAR < <(
    awk -F'\t' -v id="$TASK" '
        NR > 1 && $1 == id {
            print $2, $3, $4, $5
        }
    ' results/genotype_export_tasks.tsv
)

if [[ -z "${LANDSCAPE:-}" ]]; then
    echo "ERROR: Could not find task $TASK"
    exit 1
fi

VCF="validation_genotypes/${LANDSCAPE}/${MODEL}/${REP}/${LANDSCAPE}_${MODEL}_${REP}_year${YEAR}.vcf"

OUT="${OUTROOT}/raw/roh_task${TASKPAD}.txt"

if [[ ! -f "$VCF" ]]; then
    echo "ERROR: VCF not found: $VCF"
    exit 1
fi

echo "============================================"
echo "VALIDATION FROH"
echo "============================================"
echo "Task       : $TASK"
echo "Landscape  : $LANDSCAPE"
echo "Model      : $MODEL"
echo "Rep        : $REP"
echo "Year       : $YEAR"
echo "VCF        : $VCF"
echo "ROH output : $OUT"
echo

# ============================================================
# 2. BCFTOOLS ROH
#
# Match Nolen et al. as closely as possible:
#
# BCFtools/RoH
# genotype calls
# unseen genotype PL = 30
# fixed AF = 0.4
# ignore homozygous-reference genotypes
# Cy. semiargus recombination rate = 2.7594e-8 / bp
# ============================================================

module load bcftools/1.20

echo "BCFtools: $(bcftools --version | head -1)"

bcftools roh \
    -G30 \
    --AF-dflt 0.4 \
    --ignore-homref \
    -M "$ROH_RATE" \
    -Or \
    -o "$OUT" \
    "$VCF"

echo
echo -n "ROH regions written: "
grep -c '^RG' "$OUT" || true

# ============================================================
# 3. CALCULATE FROH >100 kb + MATCH SAMPLING DRAWS
# ============================================================

module load PDC
module load R/4.5.2-cpeGNU-26.03

Rscript "$REPO_ROOT/00-scripts/13_parse_validation_froh.R" "$TASK" "$ROH_TAG"

echo
echo "============================================"
echo "TASK COMPLETE"
echo "============================================"
