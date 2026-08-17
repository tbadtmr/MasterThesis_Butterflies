#!/bin/bash -l
#SBATCH -A lu2026-2-31
#SBATCH -p lu48
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 7-00:00:00
#SBATCH -J skane_summary
#SBATCH --array=1-50
#SBATCH -o logs/%x_%A_rep%a.out
#SBATCH -e logs/%x_%A_rep%a.err

# ============================================================================
# Cyaniris semiargus — Skåne lightweight summary simulations
#
# Computationally reduced version of the final production model used to obtain
# a larger replicate set for demographic and genomic summary trajectories.
#
# The biological model and final parameterisation are unchanged from the
# full-area production simulations, but the spatial domain is restricted to
# Skåne and large full-population SLiM snapshots are not written.
#
# Final thesis parameterisation:
#   DENS_MAX      = 120
#   K_EXPONENT    = 1.3
#   K_GLOBAL_CAP  = 120000
#   REGION_MODE   = skane_only
#
# Scenarios:
#   status_quo, restore_2km, restore_4km, restore_6km
#
# Clean reproducibility set:
#   50 replicate IDs per scenario.
#   Replicate n uses seed 100000 + n.
#
# During the thesis the simulations were distributed across HPC runs; the
# cleaned launcher covers the complete replicate range 1–50.
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SCRIPT="${SCRIPT_DIR}/simulation_skane_summaries.slim"
OUTROOT="${SCRIPT_DIR}/output"

# LUNARC SLiM installation can be supplied explicitly if needed.
SLIM_BIN="${SLIM_BIN:-$(command -v slim || true)}"

if [[ -z "$SLIM_BIN" || ! -x "$SLIM_BIN" ]]; then
    echo "ERROR: SLiM executable not found."
    echo "Set SLIM_BIN=/path/to/slim before submission if it is not in PATH."
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

if [[ -z "${BURNIN_SNAPSHOT:-}" ]]; then
    echo "ERROR: BURNIN_SNAPSHOT was not set."
    exit 1
fi

if [[ ! -f "$BURNIN_SNAPSHOT" ]]; then
    echo "ERROR: Burn-in snapshot does not exist:"
    echo "$BURNIN_SNAPSHOT"
    exit 1
fi

if [[ -z "${BURNIN_M2_FILE:-}" ]]; then
    echo "ERROR: BURNIN_M2_FILE was not set."
    exit 1
fi

if [[ ! -f "$BURNIN_M2_FILE" ]]; then
    echo "ERROR: Burn-in m2 file does not exist:"
    echo "$BURNIN_M2_FILE"
    exit 1
fi

BURNIN_SNAPSHOT="$(readlink -f "$BURNIN_SNAPSHOT")"
BURNIN_M2_FILE="$(readlink -f "$BURNIN_M2_FILE")"

# ---------------------------------------------------------------------------
# Final thesis configuration
# ---------------------------------------------------------------------------

REGION_MODE="skane_only"

DENS_MAX="120.0"
K_EXPONENT="1.3"
K_GLOBAL_CAP="120000.0"

REP=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")

# Same replicate number = same seed across scenarios.
SEED=$((100000 + SLURM_ARRAY_TASK_ID))

OUTDIR="${OUTROOT}/${SCENARIO}/rep${REP}/"

echo "========================================"
echo "C. semiargus Skåne summary simulation"
echo "Scenario       : $SCENARIO"
echo "Region         : $REGION_MODE"
echo "DENS_MAX       : $DENS_MAX"
echo "K_EXPONENT     : $K_EXPONENT"
echo "K_GLOBAL_CAP   : $K_GLOBAL_CAP"
echo "Replicate      : rep${REP}"
echo "Seed           : $SEED"
echo "Burn-in state  : $BURNIN_SNAPSHOT"
echo "Script         : $SCRIPT"
echo "Output folder  : $OUTDIR"
echo "========================================"

if [[ -d "$OUTDIR" ]] && \
   [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "ERROR: Output directory already exists and is not empty:"
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
    -d BURNIN_SNAPSHOT="\"${BURNIN_SNAPSHOT}\"" \
    -d BURNIN_M2_FILE="\"${BURNIN_M2_FILE}\"" \
    "$SCRIPT"

echo "DONE"
echo "Output: $OUTDIR"
