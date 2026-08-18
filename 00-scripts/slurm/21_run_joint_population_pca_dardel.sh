#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J popPCA
#SBATCH -o logs/popPCA_%A_%a.out
#SBATCH -e logs/popPCA_%A_%a.err

set -euo pipefail

module load bioinfo-tools
module load plink/1.90b4.9

# Required at submission.
: "${WORK_DIR:?Set WORK_DIR before submitting}"
: "${STRUCTURE_INVENTORY:?Set STRUCTURE_INVENTORY before submitting}"

WORK="$(realpath "$WORK_DIR")"
INV="$(realpath "$STRUCTURE_INVENTORY")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILDER="${JOINT_VCF_BUILDER:-$SCRIPT_DIR/../21_build_joint_population_vcf.py}"

OUT="${STRUCTURE_OUT:-$WORK/combined/results/structure}"

mkdir -p \
    "$OUT/joint_vcf" \
    "$OUT/plink" \
    "$OUT/pca"

TASK="${SLURM_ARRAY_TASK_ID}"
LINE=$(sed -n "$((TASK + 1))p" "$INV")

if [[ -z "$LINE" ]]; then
    echo "ERROR: no inventory row for array task $TASK"
    exit 1
fi

IFS=$'\t' read -r task_id rep <<< "$LINE"

echo "======================================"
echo "Replicate: $rep"
echo "======================================"

# Standard VCF layout produced by the population VCF-export workflow.
v1900="$WORK/genotypes/vcf/historical/${rep}_year1900.vcf"
v2020="$WORK/genotypes/vcf/historical/${rep}_year2020.vcf"
vsq="$WORK/genotypes/vcf/status_quo/${rep}_year2140.vcf"
vr2="$WORK/genotypes/vcf/restore_2km/${rep}_year2140.vcf"
vr4="$WORK/genotypes/vcf/restore_4km/${rep}_year2140.vcf"
vr6="$WORK/genotypes/vcf/restore_6km/${rep}_year2140.vcf"

for f in "$v1900" "$v2020" "$vsq" "$vr2" "$vr4" "$vr6"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing VCF: $f"
        exit 1
    fi
done

JOINT="$OUT/joint_vcf/${rep}_six_states.vcf"

RAW="$OUT/plink/${rep}_six_states"
PRUNE="$OUT/plink/${rep}_six_states_prune"
PRUNED="$OUT/plink/${rep}_six_states_pruned"

PCA="$OUT/pca/${rep}_six_states_PCA"

# ------------------------------------------------------------
# Build one neutral VCF shared across all six states.
# This places all states of a replicate on the same PCA axes.
# ------------------------------------------------------------

if [[ ! -s "$JOINT" ]]; then
    python3 "$BUILDER" \
        --output "$JOINT" \
        --input "y1900=$v1900" \
        --input "y2020=$v2020" \
        --input "sq2140=$vsq" \
        --input "r2_2140=$vr2" \
        --input "r4_2140=$vr4" \
        --input "r6_2140=$vr6"
fi

echo "Joint samples:"
grep -m1 '^#CHROM' "$JOINT" |
    awk -F'\t' '{print NF-9}'

echo "Joint neutral variants:"
grep -vc '^#' "$JOINT"

# ------------------------------------------------------------
# PLINK conversion
# ------------------------------------------------------------

plink \
    --vcf "$JOINT" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids '@:#' \
    --make-bed \
    --out "$RAW"

# ------------------------------------------------------------
# LD pruning: window 50 variants, step 5, r2 = 0.2
# ------------------------------------------------------------

plink \
    --bfile "$RAW" \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.2 \
    --out "$PRUNE"

NPRUNE=$(wc -l < "${PRUNE}.prune.in")

echo "Variants retained after LD pruning: $NPRUNE"

if [[ "$NPRUNE" -lt 10 ]]; then
    echo "ERROR: too few variants after LD pruning"
    exit 1
fi

plink \
    --bfile "$RAW" \
    --allow-extra-chr \
    --extract "${PRUNE}.prune.in" \
    --make-bed \
    --out "$PRUNED"

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

plink \
    --bfile "$PRUNED" \
    --allow-extra-chr \
    --pca 10 \
    --out "$PCA"

echo "PCA individuals:"
wc -l "${PCA}.eigenvec"

echo "Eigenvalues:"
head "${PCA}.eigenval"

echo "PASS"
