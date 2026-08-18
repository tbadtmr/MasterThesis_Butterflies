#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J popVCF
#SBATCH -o logs/popVCF_%A_%a.out
#SBATCH -e logs/popVCF_%A_%a.err

set -euo pipefail

# ============================================================
# Required at submission:
#
#   WORK_DIR   analysis working directory
#   TASKS_FILE task inventory with columns:
#              task_id scenario source_scenario replicate
#              year cluster snapshot
#
# Array size is supplied with sbatch --array.
# ============================================================

: "${WORK_DIR:?Set WORK_DIR before submitting}"
: "${TASKS_FILE:?Set TASKS_FILE before submitting}"

WORK="$(realpath "$WORK_DIR")"
TASKS="$(realpath "$TASKS_FILE")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Cluster-specific defaults; may be overridden at submission.
SLIM_BIN="${SLIM_BIN:-$HOME/tabea_work/conda_envs/slim5/bin/slim}"
COORD_FILE="${COORD_FILE:-$HOME/tabea_work/03_model_input/coord_chr9.txt}"
VCF_SCRIPT="${VCF_SCRIPT:-$SCRIPT_DIR/../20_export_population_manifest_vcf.slim}"

if [[ ! -x "$SLIM_BIN" ]]; then
    echo "ERROR: SLiM executable missing: $SLIM_BIN"
    exit 1
fi

for f in "$TASKS" "$COORD_FILE" "$VCF_SCRIPT"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file missing: $f"
        exit 1
    fi
done

TASK_ID="${SLURM_ARRAY_TASK_ID}"
LINE=$(sed -n "$((TASK_ID + 1))p" "$TASKS")

if [[ -z "$LINE" ]]; then
    echo "ERROR: no task row for array task $TASK_ID"
    exit 1
fi

IFS=$'\t' read -r \
    task_id scenario source_scenario replicate year cluster snapshot \
    <<< "$LINE"

MANIFEST="$WORK/sampling/manifests/$scenario/${replicate}_year${year}_manifest.tsv"

if [[ ! -f "$snapshot" ]]; then
    echo "ERROR: snapshot missing: $snapshot"
    exit 1
fi

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: manifest missing: $MANIFEST"
    exit 1
fi

OUTDIR="$WORK/genotypes/vcf/$scenario"
IDXDIR="$WORK/genotypes/indices/$scenario"

mkdir -p "$OUTDIR" "$IDXDIR"

OUT="$OUTDIR/${replicate}_year${year}.vcf"
INDICES="$IDXDIR/${replicate}_year${year}_indices.txt"

# ind_index is column 14 in the population sampling manifest.
awk -F'\t' 'NR>1 {print $14}' "$MANIFEST" > "$INDICES"

REQUEST_N=$(wc -l < "$INDICES")

if [[ "$REQUEST_N" -eq 0 ]]; then
    echo "ERROR: zero requested individuals"
    exit 1
fi

if [[ -s "$OUT" ]]; then
    EXIST_N=$(grep -m1 '^#CHROM' "$OUT" | awk -F'\t' '{print NF-9}')

    if [[ "$EXIST_N" -eq "$REQUEST_N" ]]; then
        echo "SKIP: valid VCF already exists: $OUT"
        exit 0
    fi
fi

TMP="${OUT}.tmp.${SLURM_JOB_ID}"

rm -f "$TMP"

echo "Scenario    : $scenario"
echo "Replicate   : $replicate"
echo "Year        : $year"
echo "Cluster     : $cluster"
echo "Snapshot    : $snapshot"
echo "Manifest    : $MANIFEST"
echo "Individuals : $REQUEST_N"
echo "Output      : $OUT"

"$SLIM_BIN" \
    -d SNAPSHOT="'${snapshot}'" \
    -d INDICES_FILE="'${INDICES}'" \
    -d OUT_VCF="'${TMP}'" \
    -d COORD_FILE="'${COORD_FILE}'" \
    "$VCF_SCRIPT"

VCF_N=$(grep -m1 '^#CHROM' "$TMP" | awk -F'\t' '{print NF-9}')

if [[ "$VCF_N" -ne "$REQUEST_N" ]]; then
    echo "ERROR: requested $REQUEST_N individuals but VCF contains $VCF_N"
    rm -f "$TMP"
    exit 1
fi

mv "$TMP" "$OUT"

NVAR=$(grep -vc '^#' "$OUT")

echo "Requested individuals : $REQUEST_N"
echo "VCF samples           : $VCF_N"
echo "VCF variants          : $NVAR"
echo "PASS"
