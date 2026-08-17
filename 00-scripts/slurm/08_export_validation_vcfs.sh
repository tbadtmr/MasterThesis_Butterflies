#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J valvcf
#SBATCH --array=1-94%10
#SBATCH -o valvcf_%A_%a.out
#SBATCH -e valvcf_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="${REPO_ROOT}/07_model_comparison"

cd "$WORK"

SLIM_BIN="${SLIM_BIN:-$(command -v slim || true)}"

TASKS="${WORK}/results/genotype_export_tasks.tsv"
COORD_FILE="${REPO_ROOT}/03-model_input/coord_chr9.txt"

SIM_ROOT="${SIM_ROOT:-${REPO_ROOT}/06-simulation/01-full-area-snapshots/output}"

# +1 because the task file has a header
LINE=$((SLURM_ARRAY_TASK_ID + 1))

ROW=$(sed -n "${LINE}p" "$TASKS")

if [[ -z "$ROW" ]]; then
    echo "ERROR: no task row for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

IFS=$'\t' read -r \
    TASK_ID LANDSCAPE MODEL REP YEAR REQUEST_REL \
    <<< "$ROW"

REQUEST_TSV="${WORK}/${REQUEST_REL}"

SNAPDIR="$SIM_ROOT/$LANDSCAPE/$MODEL/status_quo/$REP"

mapfile -t MATCHES < <(
    find "$SNAPDIR" \
        -maxdepth 1 \
        -type f \
        -name "*fulloutput_year${YEAR}_*.slim" \
        | sort
)

if [[ "${#MATCHES[@]}" -ne 1 ]]; then
    echo "ERROR: expected exactly one snapshot."
    echo "Directory : $SNAPDIR"
    echo "Year      : $YEAR"
    echo "Matches   : ${#MATCHES[@]}"
    printf '%s\n' "${MATCHES[@]}"
    exit 1
fi

SNAPSHOT="${MATCHES[0]}"

OUTDIR="${WORK}/validation_genotypes/${LANDSCAPE}/${MODEL}/${REP}"
mkdir -p "$OUTDIR"

PREFIX="${LANDSCAPE}_${MODEL}_${REP}_year${YEAR}"

OUT_VCF="$OUTDIR/${PREFIX}.vcf"
OUT_MAP="$OUTDIR/${PREFIX}_sample_map.tsv"

echo "============================================================"
echo "VALIDATION VCF EXPORT"
echo "============================================================"
echo "Task      : $TASK_ID / 94"
echo "Landscape : $LANDSCAPE"
echo "Model     : $MODEL"
echo "Rep       : $REP"
echo "Year      : $YEAR"
echo "Snapshot  : $SNAPSHOT"
echo "Request   : $REQUEST_TSV"
echo "VCF       : $OUT_VCF"
echo "Map       : $OUT_MAP"
echo

[[ -s "$SNAPSHOT" ]] || {
    echo "ERROR: snapshot not found or empty"
    exit 1
}

[[ -s "$REQUEST_TSV" ]] || {
    echo "ERROR: request file not found or empty"
    exit 1
}

REQUEST_N=$(( $(wc -l < "$REQUEST_TSV") - 1 ))

# Don't silently overwrite a completed valid result.
if [[ -s "$OUT_VCF" && -s "$OUT_MAP" ]]; then

    EXIST_VCF_N=$(grep -m1 '^#CHROM' "$OUT_VCF" | awk '{print NF - 9}')
    EXIST_MAP_N=$(( $(wc -l < "$OUT_MAP") - 1 ))

    if [[ "$EXIST_VCF_N" -eq "$REQUEST_N" &&
          "$EXIST_MAP_N" -eq "$REQUEST_N" ]]; then

        echo "SKIP: valid output already exists."
        echo "Request individuals : $REQUEST_N"
        echo "VCF samples         : $EXIST_VCF_N"
        echo "Map individuals     : $EXIST_MAP_N"
        exit 0
    fi
fi

# Write temporary outputs first.
TMP_VCF="${OUT_VCF}.tmp.${SLURM_JOB_ID}"
TMP_MAP="${OUT_MAP}.tmp.${SLURM_JOB_ID}"

rm -f "$TMP_VCF" "$TMP_MAP"

"$SLIM_BIN" \
    -d "SNAPSHOT='$SNAPSHOT'" \
    -d "REQUEST_TSV='$REQUEST_TSV'" \
    -d "OUT_VCF='$TMP_VCF'" \
    -d "OUT_MAP='$TMP_MAP'" \
    -d "COORD_FILE='$COORD_FILE'" \
    "$REPO_ROOT/00-scripts/08_vcf_from_validation_request.slim"

VCF_N=$(grep -m1 '^#CHROM' "$TMP_VCF" | awk '{print NF - 9}')
MAP_N=$(( $(wc -l < "$TMP_MAP") - 1 ))
VARIANT_N=$(grep -vc '^#' "$TMP_VCF")

echo
echo "QC:"
echo "Request individuals : $REQUEST_N"
echo "VCF samples         : $VCF_N"
echo "Map individuals     : $MAP_N"
echo "VCF variant records : $VARIANT_N"

if [[ "$VCF_N" -ne "$REQUEST_N" ]]; then
    echo "ERROR: VCF sample count mismatch."
    exit 1
fi

if [[ "$MAP_N" -ne "$REQUEST_N" ]]; then
    echo "ERROR: sample-map count mismatch."
    exit 1
fi

mv "$TMP_VCF" "$OUT_VCF"
mv "$TMP_MAP" "$OUT_MAP"

echo
echo "PASS"
echo "Final VCF: $OUT_VCF"
