#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Per-chromosome QC table for Cyaniris semiargus reference genome
#
# Input:
#   01-data/ncbi_dataset/data/GCA_905187585.1/
#   GCA_905187585.1_ilCyaSemi1.1_genomic.fna
#
# Output:
#   summaries/per_chr_qc.tsv
#
# Notes:
# - QC values are calculated directly from the FASTA.
# - Chromosome labels are assigned using the accession names
#   of assembly GCA_905187585.1 / ilCyaSemi1.1.
# - MT and small unplaced scaffolds are excluded from this table.
# ============================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

FA="$REPO_ROOT/01-data/ncbi_dataset/data/GCA_905187585.1/GCA_905187585.1_ilCyaSemi1.1_genomic.fna"
FAI="${FA}.fai"

OUTDIR="$REPO_ROOT/summaries"
OUT="$OUTDIR/per_chr_qc.tsv"

mkdir -p "$OUTDIR"

if [[ ! -f "$FA" ]]; then
    echo "ERROR: FASTA not found: $FA" >&2
    exit 1
fi

if [[ ! -f "$FAI" ]]; then
    echo "FASTA index not found. Creating it with samtools faidx..."
    samtools faidx "$FA"
fi

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

BASIC="$TMPDIR/basic.tsv"
NS="$TMPDIR/Ns.tsv"
DUSTMASKED="$TMPDIR/genome.dustmasked.fna"
DUST="$TMPDIR/dust.tsv"
MAXDUST="$TMPDIR/maxDustWindow_100kb.tsv"
MAP="$TMPDIR/chrom_map.tsv"

# ------------------------------------------------------------
# Chromosome label / order map for GCA_905187585.1 / ilCyaSemi1.1
# ------------------------------------------------------------
cat > "$MAP" << 'MAPEOF'
LR994546.1      1       1
LR994547.1      Z       2
LR994548.1      2       3
LR994549.1      3       4
LR994550.1      4       5
LR994551.1      5       6
LR994552.1      6       7
LR994553.1      7       8
LR994554.1      8       9
LR994555.1      9       10
LR994556.1      10      11
LR994557.1      11      12
LR994558.1      12      13
LR994559.1      13      14
LR994560.1      14      15
LR994561.1      15      16
LR994562.1      16      17
LR994563.1      17      18
LR994564.1      18      19
LR994565.1      19      20
LR994566.1      20      21
LR994567.1      21      22
LR994568.1      22      23
LR994569.1      23      24
MAPEOF

# ------------------------------------------------------------
# 1) Basic stats: accession, length, GC
# ------------------------------------------------------------
seqkit fx2tab -n -i -l -g "$FA" \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' \
> "$BASIC"

# ------------------------------------------------------------
# 2) N content per accession
# ------------------------------------------------------------
seqkit fx2tab -i -s "$FA" \
| awk 'BEGIN{OFS="\t"}{
    id=$1
    seq=$2
    len=length(seq)
    if (len==0) next

    seq_copy=seq
    n=gsub(/[Nn]/, "", seq_copy)

    print id, n, n/len
}' > "$NS"

# ------------------------------------------------------------
# 3) Low-complexity content using dustmasker
# ------------------------------------------------------------
dustmasker -in "$FA" -outfmt fasta -out "$DUSTMASKED"

seqkit fx2tab -i -s "$DUSTMASKED" \
| awk 'BEGIN{OFS="\t"}{
    id=$1
    seq=$2
    len=length(seq)
    if (len==0) next

    seq_copy=seq
    lc=gsub(/[acgtn]/, "", seq_copy)

    print id, lc, lc/len
}' > "$DUST"

# ------------------------------------------------------------
# 4) Maximum low-complexity fraction in any 100 kb window
# ------------------------------------------------------------
seqkit sliding -s 100000 -W 100000 "$DUSTMASKED" \
| seqkit fx2tab -i -s \
| awk 'BEGIN{OFS="\t"}{
    id=$1
    seq=$2
    len=length(seq)
    if (len==0) next

    split(id, a, ":")
    acc=a[1]
    sub(/_sliding$/, "", acc)

    seq_copy=seq
    lc=gsub(/[acgtn]/, "", seq_copy)
    frac=lc/len

    if (!(acc in max) || frac > max[acc]) {
        max[acc] = frac
    }
}
END{
    for (acc in max) print acc, max[acc]
}' > "$MAXDUST"

# ------------------------------------------------------------
# 5) Merge all tables and keep only mapped chromosomes
# ------------------------------------------------------------
awk 'BEGIN{OFS="\t"}
FNR==NR {
    chr[$1] = $2
    ord[$1] = $3
    keep[$1] = 1
    next
}
FILENAME==ARGV[2] {
    len[$1] = $2
    gc[$1]  = $3
    next
}
FILENAME==ARGV[3] {
    ncount[$1] = $2
    nfrac[$1]  = $3
    next
}
FILENAME==ARGV[4] {
    lcbases[$1] = $2
    lcfrac[$1]  = $3
    next
}
FILENAME==ARGV[5] {
    maxlc[$1] = $2
    next
}
END {
    print "Chr_number","Accession","Size_bp","GC_percent","N_count","N_fraction","LowComplexity_bases","LowComplexity_fraction","MaxLowComplexity_100kb"

    for (acc in keep) {
        print chr[acc], acc, len[acc], gc[acc], ncount[acc], nfrac[acc], lcbases[acc], lcfrac[acc], maxlc[acc]
    }
}' "$MAP" "$BASIC" "$NS" "$DUST" "$MAXDUST" \
| { read header; echo "$header"; sort -k1,1n; } \
> "$OUT"

echo "QC table written to: $OUT"
echo
column -t "$OUT" | head -15
