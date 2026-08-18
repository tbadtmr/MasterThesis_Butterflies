#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 04:00:00
#SBATCH -J skADM
#SBATCH --array=1-144
#SBATCH -o combined/results/structure/logs/skADM_%A_%a.out
#SBATCH -e combined/results/structure/logs/skADM_%A_%a.err

set -euo pipefail

module load bioinfo-tools
module load ADMIXTURE/1.3.0

WORK="$HOME/tabea_work/08_population_analysis_final/skane_only"

TASKS="$WORK/combined/inventories/admixture_tasks.tsv"
PLINKROOT="$WORK/combined/results/structure/plink"
OUTROOT="$WORK/combined/results/structure/admixture"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TASKS")

IFS=$'\t' read -r task_id REP K <<< "$LINE"

PREFIX="$PLINKROOT/${REP}_six_states_pruned"

for ext in bed bim fam; do
    if [[ ! -s "${PREFIX}.${ext}" ]]; then
        echo "ERROR: missing ${PREFIX}.${ext}"
        exit 1
    fi
done

OUTDIR="$OUTROOT/$REP/K${K}"
mkdir -p "$OUTDIR"

REP_NUM=$((10#${REP#rep}))

echo "=========================================="
echo "Replicate : $REP"
echo "K         : $K"
echo "Starts    : 20"
echo "CV folds  : 5"
echo "=========================================="

# ------------------------------------------------------------
# 20 independent starts
# ------------------------------------------------------------

for RUN in $(seq 1 20); do

    RUNDIR=$(printf "%s/run%02d" "$OUTDIR" "$RUN")
    mkdir -p "$RUNDIR"

    cd "$RUNDIR"

    # ADMIXTURE expects the .bed/.bim/.fam files to share a basename.
    ln -sf "${PREFIX}.bed" data.bed
    ln -sf "${PREFIX}.fam" data.fam

    # ADMIXTURE requires integer chromosome codes.
    # All simulated variants are from chromosome 9.
    awk 'BEGIN{OFS="\t"} {$1=9; print}'         "${PREFIX}.bim" > data.bim

    # Fixed reproducible but different seed for every
    # replicate × K × start combination.
    SEED=$((20260817 + REP_NUM * 10000 + K * 100 + RUN))

    echo
    echo "RUN $RUN / seed $SEED"

    admixture \
        --cv=5 \
        -j4 \
        -s "$SEED" \
        data.bed \
        "$K" \
        > run.log 2>&1

    if [[ ! -s "data.${K}.Q" ]]; then
        echo "ERROR: missing Q output for run $RUN"
        exit 1
    fi

    if [[ ! -s "data.${K}.P" ]]; then
        echo "ERROR: missing P output for run $RUN"
        exit 1
    fi

done

echo
echo "=========================================="
echo "ALL 20 RUNS COMPLETE"
echo "Replicate: $REP"
echo "K: $K"
echo "=========================================="
