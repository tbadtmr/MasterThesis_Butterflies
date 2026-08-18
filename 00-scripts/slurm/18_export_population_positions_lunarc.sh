#!/bin/bash -l
#SBATCH -A lu2026-2-31
#SBATCH -p lu48
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J popPOS
#SBATCH -o logs/popPOS_%A_%a.out
#SBATCH -e logs/popPOS_%A_%a.err

set -euo pipefail

: "${WORK_DIR:?Set WORK_DIR before submitting}"
: "${TASKS_FILE:?Set TASKS_FILE before submitting}"

WORK="$(realpath "$WORK_DIR")"
TASKS="$(realpath "$TASKS_FILE")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BASE="${LUNARC_BASE:-$HOME/lu2026-12-12/Tabea}"

SLIM_BIN="${SLIM_BIN:-$BASE/tools/mamba_root/envs/slim5/bin/slim}"
COORD_FILE="${COORD_FILE:-$BASE/03_model_input/coord_chr9.txt}"
POSITION_SCRIPT="${POSITION_SCRIPT:-$SCRIPT_DIR/../18_export_population_positions.slim}"

for f in "$TASKS" "$COORD_FILE" "$POSITION_SCRIPT"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file missing: $f"
        exit 1
    fi
done

if [[ ! -x "$SLIM_BIN" ]]; then
    echo "ERROR: SLiM executable missing: $SLIM_BIN"
    exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID}"
LINE=$(sed -n "$((TASK_ID + 1))p" "$TASKS")

if [[ -z "$LINE" ]]; then
    echo "ERROR: no task row for array task $TASK_ID"
    exit 1
fi

IFS=$'\t' read -r \
    task_id scenario source_scenario replicate year cluster snapshot \
    <<< "$LINE"

if [[ ! -f "$snapshot" ]]; then
    echo "ERROR: snapshot missing: $snapshot"
    exit 1
fi

OUTDIR="$WORK/positions/$scenario"
mkdir -p "$OUTDIR"

OUT="$OUTDIR/${replicate}_year${year}_positions.tsv"

if [[ -s "$OUT" ]]; then
    echo "SKIP: position file already exists: $OUT"
    exit 0
fi

TMP="${OUT}.tmp.${SLURM_JOB_ID}"

rm -f "$TMP"

echo "Scenario   : $scenario"
echo "Source     : $source_scenario"
echo "Replicate  : $replicate"
echo "Year       : $year"
echo "Cluster    : $cluster"
echo "Snapshot   : $snapshot"
echo "Output     : $OUT"

"$SLIM_BIN" \
    -d SNAPSHOT="'${snapshot}'" \
    -d OUT="'${TMP}'" \
    -d MODEL="'${scenario}'" \
    -d REP="'${replicate}'" \
    -d YEAR="'${year}'" \
    -d COORD_FILE="'${COORD_FILE}'" \
    "$POSITION_SCRIPT"

if [[ ! -s "$TMP" ]]; then
    echo "ERROR: position export produced no output"
    rm -f "$TMP"
    exit 1
fi

mv "$TMP" "$OUT"

N=$(awk 'END{print NR-1}' "$OUT")

echo "Individuals exported: $N"
echo "PASS"
