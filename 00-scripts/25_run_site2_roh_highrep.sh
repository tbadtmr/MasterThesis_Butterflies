#!/bin/bash -l
#SBATCH -A naiss2026-1-7
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J rohSite2
#SBATCH -o logs/regional_roh_site2/roh_%A_%a.out
#SBATCH -e logs/regional_roh_site2/roh_%A_%a.err

set -euo pipefail

BASE="$HOME/tabea_work/08_population_analysis_final/skane_only"
TASKS="$BASE/inventories/roh_site2_tasks.tsv"
OUTROOT="$BASE/results/regional_roh_site2"

PARSER="$HOME/tabea_work/08_population_analysis/scripts/26_parse_site2_roh_highrep.R"

mkdir -p \
    "$OUTROOT/raw" \
    "$OUTROOT/subset_vcf" \
    "$OUTROOT/manifests" \
    "$OUTROOT/segments" \
    "$OUTROOT/individuals" \
    "$OUTROOT/summaries"

LINE=$(awk -F'\t' \
    -v n="$SLURM_ARRAY_TASK_ID" \
    'NR == n + 1 {print; exit}' \
    "$TASKS")

if [[ -z "$LINE" ]]; then
    echo "ERROR: no task for array index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

IFS=$'\t' read -r \
    SCENARIO REP YEAR VCF MANIFEST \
    <<< "$LINE"

BASE_NAME="${SCENARIO}_${REP}_year${YEAR}_site2"

MANIFEST_DIR="$OUTROOT/manifests/$SCENARIO"
RAW_DIR="$OUTROOT/raw/$SCENARIO"
VCF_OUT_DIR="$OUTROOT/subset_vcf/$SCENARIO"

mkdir -p \
    "$MANIFEST_DIR" \
    "$RAW_DIR" \
    "$VCF_OUT_DIR"

SITE_MANIFEST="$MANIFEST_DIR/${BASE_NAME}_manifest.tsv"
SAMPLES="$MANIFEST_DIR/${BASE_NAME}_samples.txt"

SUBVCF="$VCF_OUT_DIR/${BASE_NAME}.vcf"
RAW="$RAW_DIR/${BASE_NAME}_roh.txt"

# ============================================================
# 1. Keep only fixed sampling site number 2
# ============================================================

awk -F'\t' '
BEGIN { OFS="\t" }

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "site_number") sn = i
    }

    if (!sn) {
        print "ERROR: site_number column missing" > "/dev/stderr"
        exit 2
    }

    print
    next
}

$sn == 2 {
    print
}
' "$MANIFEST" > "$SITE_MANIFEST"

N_MAN=$(( $(wc -l < "$SITE_MANIFEST") - 1 ))

if [[ "$N_MAN" -lt 1 ]]; then
    echo "ERROR: empty site-2 manifest" >&2
    exit 1
fi

# ============================================================
# 2. Exact VCF sample names
# ============================================================

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "ind_index") ii = i
    }

    if (!ii) {
        print "ERROR: ind_index column missing" > "/dev/stderr"
        exit 2
    }

    next
}

{
    print "p1:i" $ii
}
' "$SITE_MANIFEST" > "$SAMPLES"

# ============================================================
# 3. Subset VCF columns without rewriting VCF annotations
#
# The SLiM VCF contains header tags that bcftools can read but
# cannot safely rewrite during sample subsetting. We therefore
# retain the original VCF text and select only the requested
# sample columns.
# ============================================================

python3 - "$VCF" "$SAMPLES" "$SUBVCF" <<'PY'
import sys

vcf, sample_file, out = sys.argv[1:]

with open(sample_file) as f:
    wanted = [x.strip() for x in f if x.strip()]

wanted_set = set(wanted)

with open(vcf) as inp, open(out, "w") as dst:

    sample_indices = None
    selected_names = None

    for line in inp:

        if line.startswith("##"):
            dst.write(line)
            continue

        if line.startswith("#CHROM"):

            fields = line.rstrip("\n").split("\t")
            samples = fields[9:]

            lookup = {
                sample: 9 + i
                for i, sample in enumerate(samples)
            }

            missing = [
                s for s in wanted
                if s not in lookup
            ]

            if missing:
                raise SystemExit(
                    "Samples missing from VCF: "
                    + ", ".join(missing[:20])
                )

            # Preserve requested manifest order.
            sample_indices = [
                lookup[s]
                for s in wanted
            ]

            selected_names = [
                fields[i]
                for i in sample_indices
            ]

            dst.write(
                "\t".join(
                    fields[:9] + selected_names
                )
                + "\n"
            )

            continue

        if line.startswith("#"):
            dst.write(line)
            continue

        if sample_indices is None:
            raise SystemExit(
                "No #CHROM header found before variant records."
            )

        fields = line.rstrip("\n").split("\t")

        dst.write(
            "\t".join(
                fields[:9] +
                [fields[i] for i in sample_indices]
            )
            + "\n"
        )

if selected_names is None:
    raise SystemExit("No #CHROM line found.")

if len(selected_names) != len(wanted):
    raise SystemExit(
        "Subset sample count does not match requested count."
    )

print(
    f"Subset VCF: {len(selected_names)} samples"
)
PY

N_VCF=$(awk -F'\t' \
    '/^#CHROM/ {print NF-9; exit}' \
    "$SUBVCF")

if [[ "$N_VCF" -ne "$N_MAN" ]]; then
    echo "ERROR: VCF/manifest sample mismatch: $N_VCF vs $N_MAN" >&2
    exit 1
fi

echo "=============================================="
echo "SITE-2 REGIONAL ROH"
echo "=============================================="
echo "Scenario    : $SCENARIO"
echo "Replicate   : $REP"
echo "Year        : $YEAR"
echo "Individuals : $N_VCF"
echo "=============================================="

# ============================================================
# 4. BCFtools/RoH
#
# Same settings as validated ROH analysis
# ============================================================

module load bcftools/1.20

bcftools roh \
    -G30 \
    --AF-dflt 0.4 \
    --ignore-homref \
    -M 2.7594e-8 \
    -Or \
    -o "$RAW" \
    "$SUBVCF"

# ============================================================
# 5. Parse >100-kb ROH
# ============================================================

module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

Rscript "$PARSER" \
    "$RAW" \
    "$SITE_MANIFEST" \
    "$OUTROOT" \
    "$BASE_NAME"

echo "DONE: $BASE_NAME"
