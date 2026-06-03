#!/usr/bin/env bash
set -euo pipefail

# Make continuous exon/intron/non-coding coordination file for chromosome 9.
# Output coordinates are 1-based inclusive.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ASM_DIR="$REPO_ROOT/01-data/ncbi_dataset/data/GCA_905187585.1"
FASTA="$ASM_DIR/GCA_905187585.1_ilCyaSemi1.1_genomic.fna"
FAI="${FASTA}.fai"
GFF3="$ASM_DIR/genes.gff3"

OUTDIR="$REPO_ROOT/03-model_input"
TESTDIR="$OUTDIR/testing"

CHR_LABEL="9"
CHR_ACC="LR994555.1"
OUTBASE="coord_chr9"
SOURCE_LABEL="braker_augustus"

mkdir -p "$OUTDIR" "$TESTDIR"

[[ -f "$FASTA" ]] || { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
[[ -f "$GFF3" ]] || { echo "ERROR: GFF3 not found: $GFF3" >&2; exit 1; }

if [[ ! -f "$FAI" ]]; then
    echo "FASTA index not found. Creating it..."
    samtools faidx "$FASTA"
fi

CHR_LEN=$(awk -v acc="$CHR_ACC" '$1==acc {print $2; exit}' "$FAI")

if [[ -z "${CHR_LEN:-}" ]]; then
    echo "ERROR: Could not find $CHR_ACC in FASTA index: $FAI" >&2
    exit 1
fi

echo "Using FASTA : $FASTA"
echo "Using GFF3  : $GFF3"
echo "Chromosome : $CHR_LABEL / $CHR_ACC"
echo "Length     : $CHR_LEN bp"

CHR_GFF="$OUTDIR/${OUTBASE}.gff3"

awk -F'\t' -v chr="$CHR_LABEL" '
    $0 !~ /^#/ && $1 == chr
' "$GFF3" > "$CHR_GFF"

if [[ ! -s "$CHR_GFF" ]]; then
    echo "ERROR: No annotation rows found for chromosome $CHR_LABEL in $GFF3" >&2
    exit 1
fi

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

GENES_BED="$TMPDIR/genes.bed"
EXONS_BED="$TMPDIR/exons.bed"
EXONS_MERGED="$TMPDIR/exons_merged.bed"
BREAKS_SORTED="$TMPDIR/breakpoints_sorted.txt"

# Internally use 0-based half-open intervals for clean interval arithmetic.
awk -F'\t' '
BEGIN { OFS="\t" }
$3 == "gene" {
    start0 = $4 - 1
    end = $5
    if (start0 < 0) start0 = 0
    print $1, start0, end
}
' "$CHR_GFF" | sort -k2,2n -k3,3n > "$GENES_BED"

awk -F'\t' '
BEGIN { OFS="\t" }
$3 == "exon" {
    start0 = $4 - 1
    end = $5
    if (start0 < 0) start0 = 0
    print $1, start0, end
}
' "$CHR_GFF" | sort -k2,2n -k3,3n > "$EXONS_BED"

if [[ ! -s "$EXONS_BED" ]]; then
    echo "ERROR: No exon rows found for chromosome $CHR_LABEL" >&2
    exit 1
fi

awk '
BEGIN { OFS="\t" }
NR == 1 {
    chr=$1; s=$2; e=$3
    next
}
{
    if ($2 <= e) {
        if ($3 > e) e=$3
    } else {
        print chr, s, e
        chr=$1; s=$2; e=$3
    }
}
END {
    if (NR > 0) print chr, s, e
}
' "$EXONS_BED" > "$EXONS_MERGED"

{
    echo 0
    echo "$CHR_LEN"
    awk '{print $2; print $3}' "$GENES_BED"
    awk '{print $2; print $3}' "$EXONS_MERGED"
} | sort -n | uniq > "$BREAKS_SORTED"

awk -v chr_label="$CHR_LABEL" -v source_label="$SOURCE_LABEL" '
BEGIN { OFS="\t" }

FNR==NR {
    gstart[++ng] = $2
    gend[ng] = $3
    next
}

ARGIND==2 {
    estart[++ne] = $2
    eend[ne] = $3
    next
}

ARGIND==3 {
    pts[++n] = $1
    next
}

END {
    seg_n = 0

    for (i=1; i<n; i++) {
        s = pts[i]
        e = pts[i+1]
        if (e <= s) continue

        in_exon = 0
        for (j=1; j<=ne; j++) {
            if (s >= estart[j] && e <= eend[j]) {
                in_exon = 1
                break
            }
        }

        in_gene = 0
        for (j=1; j<=ng; j++) {
            if (s >= gstart[j] && e <= gend[j]) {
                in_gene = 1
                break
            }
        }

        if (in_exon) {
            type = "exon"
        } else if (in_gene) {
            type = "intron"
        } else {
            type = "non_coding"
        }

        if (seg_n == 0 || s != seg_e[seg_n] || type != seg_t[seg_n]) {
            seg_n++
            seg_s[seg_n] = s
            seg_e[seg_n] = e
            seg_t[seg_n] = type
        } else {
            seg_e[seg_n] = e
        }
    }

    print "V1","V2.x","V3","V4","V5","grp","DUP","len"

    for (k=1; k<=seg_n; k++) {
        # Convert from internal 0-based half-open to output 1-based inclusive.
        start_out = seg_s[k] + 1
        end_out = seg_e[k]
        len = end_out - start_out + 1

        print chr_label, source_label, seg_t[k], start_out, end_out, k, "NO", len
    }
}
' "$GENES_BED" "$EXONS_MERGED" "$BREAKS_SORTED" > "$OUTDIR/${OUTBASE}.txt"

MID=$(( CHR_LEN / 2 ))

MID200_START=$(( MID - 100000 + 1 ))
MID200_END=$(( MID + 100000 ))

MID2M_START=$(( MID - 1000000 + 1 ))
MID2M_END=$(( MID + 1000000 ))

awk -v s="$MID200_START" -v e="$MID200_END" \
    'NR==1 || ($5>=s && $4<=e)' \
    "$OUTDIR/${OUTBASE}.txt" > "$TESTDIR/${OUTBASE}_win_mid_200kb.txt"

awk -v s="$MID2M_START" -v e="$MID2M_END" \
    'NR==1 || ($5>=s && $4<=e)' \
    "$OUTDIR/${OUTBASE}.txt" > "$TESTDIR/${OUTBASE}_win_mid_2Mb.txt"

echo
echo "Wrote:"
echo "  $OUTDIR/${OUTBASE}.gff3"
echo "  $OUTDIR/${OUTBASE}.txt"
echo "  $TESTDIR/${OUTBASE}_win_mid_200kb.txt"
echo "  $TESTDIR/${OUTBASE}_win_mid_2Mb.txt"
