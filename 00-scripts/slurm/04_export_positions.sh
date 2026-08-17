#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 06:00:00
#SBATCH -J all_positions
#SBATCH -o all_positions_%j.out
#SBATCH -e all_positions_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="${REPO_ROOT}/07_model_comparison"


# Same SLiM environment as the successful simulation jobs
SLIM_BIN="${SLIM_BIN:-$(command -v slim || true)}"

if [[ ! -x "$SLIM_BIN" ]]; then
    echo "ERROR: SLiM executable not found:"
    echo "$SLIM_BIN"
    exit 1
fi

MODEL="Kexp1p3_D120"

SIM_ROOT="${SIM_ROOT:-${REPO_ROOT}/06-simulation/01-full-area-snapshots/output}"
BASE="${SIM_ROOT}/status_quo"

# WORK defined from repository location above
SLIM_SCRIPT="${REPO_ROOT}/00-scripts/04_export_positions.slim"
COORD="${REPO_ROOT}/03-model_input/coord_chr9.txt"
OUTDIR="${WORK}/positions"

mkdir -p "$OUTDIR"

echo "=============================================="
echo "Position export"
echo "Model : $MODEL"
echo "Host  : $(hostname)"
echo "Start : $(date)"
echo "SLiM  : $SLIM_BIN"
"$SLIM_BIN" --version
echo "=============================================="

for repnum in 1 2 3 4 5; do

    REP=$(printf "rep%03d" "$repnum")

    for YEAR in 1900 1951 1956 2020; do

        SNAP=$(find "${BASE}/${REP}" \
            -maxdepth 1 \
            -type f \
            -name "*fulloutput_year${YEAR}_*.slim" \
            -print -quit)

        if [[ -z "$SNAP" ]]; then
            echo "ERROR: no snapshot for ${REP}, year ${YEAR}" >&2
            exit 1
        fi

        OUT="${OUTDIR}/${MODEL}_${REP}_year${YEAR}_positions.tsv"

        # Allows safe restart of this batch job.
        # Existing completed position files are not regenerated.
        if [[ -s "$OUT" ]]; then
            echo
            echo "SKIP: already exists"
            echo "$OUT"
            continue
        fi

        # Write to temporary file first so an interrupted export
        # cannot be mistaken for a completed output.
        TMP="${OUT}.tmp.$$"
        rm -f "$TMP"

        echo
        echo "================================================"
        echo "Model     : $MODEL"
        echo "Replicate : $REP"
        echo "Year      : $YEAR"
        echo "Snapshot  : $SNAP"
        echo "Output    : $OUT"
        echo "Started   : $(date)"
        echo "================================================"

        "$SLIM_BIN" \
            -d SNAPSHOT="'${SNAP}'" \
            -d OUT="'${TMP}'" \
            -d MODEL="'${MODEL}'" \
            -d REP="'${REP}'" \
            -d YEAR=${YEAR} \
            -d COORD_FILE="'${COORD}'" \
            "$SLIM_SCRIPT"

        mv "$TMP" "$OUT"

        echo "Completed : $(date)"
        wc -l "$OUT"

    done
done

echo
echo "=============================================="
echo "ALL POSITION EXPORTS FINISHED"
echo "Finished: $(date)"
echo "=============================================="

ls -lh "${OUTDIR}"/${MODEL}_rep*_year*_positions.tsv
