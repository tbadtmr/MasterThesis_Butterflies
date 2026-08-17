#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 7-00:00:00
#SBATCH -J full_resume
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

# ============================================================================
# Cyaniris semiargus — resume a full-area production simulation
#
# Continues one interrupted full-area replicate from a saved future snapshot.
# This was used when individual simulations exceeded the Dardel wall-time.
#
# Required at submission:
#   SCENARIO
#   REP_NUMBER
#   RESUME_SNAPSHOT
#
# RESUME_M2_FILE is optional. If omitted, the script derives the matching
# main_chr9_m2_*.txt filename from the supplied fulloutput snapshot.
#
# Final thesis parameterisation:
#   DENS_MAX      = 120
#   K_EXPONENT    = 1.3
#   K_GLOBAL_CAP  = 120000
#
# The continuation is written under output/resume/. These continuation files
# can subsequently be combined with the corresponding original replicate.
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SCRIPT="${SCRIPT_DIR}/simulation_resume_from_snapshot.slim"
OUTROOT="${SCRIPT_DIR}/output/resume"

SLIM_BIN="${SLIM_BIN:-$HOME/tabea_work/conda_envs/slim5/bin/slim}"

if [[ ! -x "$SLIM_BIN" ]]; then
    echo "ERROR: SLiM executable not found or not executable:"
    echo "$SLIM_BIN"
    exit 1
fi

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: SLiM script not found:"
    echo "$SCRIPT"
    exit 1
fi

# ---------------------------------------------------------------------------
# Required inputs
# ---------------------------------------------------------------------------

if [[ -z "${SCENARIO:-}" ]]; then
    echo "ERROR: SCENARIO was not set."
    echo "Allowed values: status_quo, restore_2km, restore_4km, restore_6km"
    exit 1
fi

case "$SCENARIO" in
    status_quo|restore_2km|restore_4km|restore_6km)
        ;;
    *)
        echo "ERROR: Unknown scenario: $SCENARIO"
        exit 1
        ;;
esac

if [[ -z "${REP_NUMBER:-}" ]]; then
    echo "ERROR: REP_NUMBER was not set."
    echo "Use an integer from 1 to 5."
    exit 1
fi

if ! [[ "$REP_NUMBER" =~ ^[1-5]$ ]]; then
    echo "ERROR: REP_NUMBER must be an integer from 1 to 5."
    exit 1
fi

if [[ -z "${RESUME_SNAPSHOT:-}" ]]; then
    echo "ERROR: RESUME_SNAPSHOT was not set."
    exit 1
fi

if [[ ! -f "$RESUME_SNAPSHOT" ]]; then
    echo "ERROR: Resume snapshot does not exist:"
    echo "$RESUME_SNAPSHOT"
    exit 1
fi

RESUME_SNAPSHOT="$(readlink -f "$RESUME_SNAPSHOT")"

# Derive the matching m2 filename unless supplied explicitly.
if [[ -z "${RESUME_M2_FILE:-}" ]]; then
    RESUME_M2_FILE="${RESUME_SNAPSHOT/fulloutput/m2}"
    RESUME_M2_FILE="${RESUME_M2_FILE%.slim}.txt"
fi

if [[ ! -f "$RESUME_M2_FILE" ]]; then
    echo "ERROR: Matching resume m2 file does not exist:"
    echo "$RESUME_M2_FILE"
    exit 1
fi

RESUME_M2_FILE="$(readlink -f "$RESUME_M2_FILE")"

# ---------------------------------------------------------------------------
# Final thesis configuration
# ---------------------------------------------------------------------------

REGION_MODE="full"

DENS_MAX="120.0"
K_EXPONENT="1.3"
K_GLOBAL_CAP="120000.0"

REP=$(printf "%03d" "$REP_NUMBER")
SEED=$((100000 + REP_NUMBER))

OUTDIR="${OUTROOT}/${SCENARIO}/rep${REP}/"

echo "========================================"
echo "C. semiargus full-area continuation"
echo "Scenario        : $SCENARIO"
echo "Region          : $REGION_MODE"
echo "Replicate       : rep${REP}"
echo "Seed            : $SEED"
echo "DENS_MAX        : $DENS_MAX"
echo "K_EXPONENT      : $K_EXPONENT"
echo "K_GLOBAL_CAP    : $K_GLOBAL_CAP"
echo "Resume snapshot : $RESUME_SNAPSHOT"
echo "Resume m2       : $RESUME_M2_FILE"
echo "Script          : $SCRIPT"
echo "Output folder   : $OUTDIR"
echo "========================================"

if [[ -d "$OUTDIR" ]] && \
   [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "ERROR: Resume output directory already exists and is not empty:"
    echo "$OUTDIR"
    exit 1
fi

mkdir -p "$OUTDIR"

"$SLIM_BIN" --version

"$SLIM_BIN" \
    -s "$SEED" \
    -d SCENARIO="\"${SCENARIO}\"" \
    -d REGION_MODE="\"${REGION_MODE}\"" \
    -d DENS_MAX="${DENS_MAX}" \
    -d K_EXPONENT="${K_EXPONENT}" \
    -d K_GLOBAL_CAP="${K_GLOBAL_CAP}" \
    -d WRITE_IBD=F \
    -d REP_ID="\"rep${REP}\"" \
    -d OUT_DIR="\"${OUTDIR}\"" \
    -d BASE_DIR="\"${REPO_ROOT}\"" \
    -d RESUME_SNAPSHOT="\"${RESUME_SNAPSHOT}\"" \
    -d RESUME_M2_FILE="\"${RESUME_M2_FILE}\"" \
    "$SCRIPT"

echo "DONE"
echo "Output: $OUTDIR"
