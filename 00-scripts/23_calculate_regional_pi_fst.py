#!/usr/bin/env python3

# ============================================================
# Regional neutral diversity and differentiation
#
# Input:
#   30 full-landscape regional VCFs
#
# Sampling:
#   7 regions
#   3 sites / region
#   target 20 individuals / site
#   = target 60 individuals / region; fewer where populations decline
#
# Metrics:
#   - neutral nucleotide diversity (pi), MT=1 only
#   - pairwise Hudson FST, MT=1 only
#
# Comparisons:
#   1900 -> 2020
#   2020 -> 2140 under each scenario
#   restoration vs paired status quo at 2140
# ============================================================

from pathlib import Path
import argparse
import csv
import math
import re
import statistics
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Calculate regional neutral pi and pairwise Hudson FST."
)

parser.add_argument(
    "--work-dir",
    required=True,
    help="Population-analysis working directory containing genotypes/ and sampling/."
)

parser.add_argument(
    "--design",
    required=True,
    help="Fixed regional sampling design TSV."
)

parser.add_argument(
    "--replicate-inventory",
    required=True,
    help="TSV containing a replicate column."
)

parser.add_argument(
    "--coord-file",
    required=True,
    help="Chromosome coordination file used by the simulation."
)

parser.add_argument(
    "--out-dir",
    default=None,
    help="Output directory (default: WORK_DIR/results/regional_genetics)."
)

args = parser.parse_args()

BASE = Path(args.work_dir).expanduser().resolve()
DESIGN_FILE = Path(args.design).expanduser().resolve()
REPLICATE_FILE = Path(args.replicate_inventory).expanduser().resolve()
COORD_FILE = Path(args.coord_file).expanduser().resolve()

VCF_DIR = BASE / "genotypes" / "vcf"
MANIFEST_DIR = BASE / "sampling" / "manifests"

OUT_DIR = (
    Path(args.out_dir).expanduser().resolve()
    if args.out_dir
    else BASE / "results" / "regional_genetics"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


def read_region_order(path):

    regions = []

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="	")

        if "region" not in (reader.fieldnames or []):
            raise RuntimeError(
                f"Design file lacks a region column: {path}"
            )

        for row in reader:
            region = row["region"]

            if region and region not in regions:
                regions.append(region)

    if not regions:
        raise RuntimeError(
            f"No regions found in design file: {path}"
        )

    return regions


def read_replicates(path):

    reps = []

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="	")

        if "replicate" not in (reader.fieldnames or []):
            raise RuntimeError(
                f"Replicate inventory lacks a replicate column: {path}"
            )

        for row in reader:
            rep = row["replicate"]

            if rep and rep not in reps:
                reps.append(rep)

    if not reps:
        raise RuntimeError(
            f"No replicates found in inventory: {path}"
        )

    return reps


REGION_ORDER = read_region_order(DESIGN_FILE)
REPLICATES = read_replicates(REPLICATE_FILE)

SCENARIO_ORDER = [
    "shared_history",
    "status_quo",
    "restore_2km",
    "restore_4km",
    "restore_6km",
]


# ============================================================
# Chromosome length
# ============================================================

def get_chrom_length(path):

    max_end = -1

    with open(path) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 5:
                continue

            try:
                end = int(fields[4])
            except ValueError:
                continue

            max_end = max(
                max_end,
                end
            )

    if max_end < 0:
        raise RuntimeError(
            f"Could not determine chromosome length from {path}"
        )

    # SLiM coordinates are 0-based inclusive
    return max_end + 1


CHROM_LEN = get_chrom_length(
    COORD_FILE
)

print(
    "Chromosome length:",
    CHROM_LEN
)


# ============================================================
# Helpers
# ============================================================

def sample_sd(values):

    if len(values) < 2:
        return float("nan")

    return statistics.stdev(
        values
    )


def mean(values):

    if not values:
        return float("nan")

    return statistics.mean(
        values
    )


def safe_ratio(a, b):

    if b == 0:
        return float("nan")

    return a / b


def parse_info(info):

    result = {}

    for item in info.split(";"):

        if "=" in item:

            key, value = item.split(
                "=",
                1
            )

            result[key] = value

        else:

            result[item] = True

    return result


def parse_gt(gt):

    # Remove any additional FORMAT information
    gt = gt.split(":")[0]

    if gt in (
        ".",
        "./.",
        ".|.",
        "~"
    ):
        return []

    gt = gt.replace(
        "|",
        "/"
    )

    alleles = []

    for a in gt.split("/"):

        if a in (
            ".",
            "~",
            ""
        ):
            continue

        alleles.append(
            int(a)
        )

    return alleles


def get_metadata_from_filename(path):

    storage_scenario = path.parent.name

    allowed = {
        "historical",
        "shared_history",
        "status_quo",
        "restore_2km",
        "restore_4km",
        "restore_6km",
    }

    if storage_scenario not in allowed:
        raise RuntimeError(
            f"Unexpected scenario directory: {storage_scenario}"
        )

    m = re.match(
        r"^(rep\\d{3})_year(\\d{4})\\.vcf$",
        path.name
    )

    if not m:
        raise RuntimeError(
            f"Unexpected VCF filename: {path.name}"
        )

    rep = m.group(1)
    year = int(m.group(2))

    # Historical VCFs are stored under historical/, while the
    # existing summary code uses shared_history internally.
    scenario = (
        "shared_history"
        if storage_scenario == "historical"
        else storage_scenario
    )

    return scenario, storage_scenario, rep, year


def read_manifest(path):

    index_meta = {}

    with open(path) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            idx = int(
                row["ind_index"]
            )

            if idx in index_meta:
                raise RuntimeError(
                    f"Duplicate ind_index {idx} in {path}"
                )

            index_meta[idx] = {
                "region": row["region"],
                "site": row["site"]
            }

    if len(index_meta) == 0:
        raise RuntimeError(
            f"Manifest contains no individuals: {path}"
        )

    invalid_regions = sorted(
        {
            x["region"]
            for x in index_meta.values()
        } - set(REGION_ORDER)
    )

    if invalid_regions:
        raise RuntimeError(
            f"Unexpected region labels in {path}: "
            f"{invalid_regions}"
        )

    return index_meta


def sample_name_to_index(name):

    # Expected:
    # p1:i14425

    m = re.search(
        r":i(\d+)$",
        name
    )

    if not m:
        raise RuntimeError(
            f"Cannot parse SLiM sample ID: {name}"
        )

    return int(
        m.group(1)
    )


# ============================================================
# Analyse one VCF
# ============================================================

def analyse_vcf(vcf_path):

    scenario, storage_scenario, rep, year = (
        get_metadata_from_filename(
            vcf_path
        )
    )

    manifest_path = (
        MANIFEST_DIR /
        storage_scenario /
        f"{rep}_year{year}_manifest.tsv"
    )

    if not manifest_path.exists():

        raise FileNotFoundError(
            f"Manifest missing for {vcf_path}\n"
            f"Expected: {manifest_path}"
        )

    meta = read_manifest(
        manifest_path
    )

    sample_regions = None
    sample_sites = None

    n_individuals_region = defaultdict(int)

    # Sum of per-site pairwise diversity
    pi_sum = defaultdict(float)

    # Hudson numerator / denominator sums
    fst_num_sum = defaultdict(float)
    fst_den_sum = defaultdict(float)
    fst_sites = defaultdict(int)

    total_variants = 0
    neutral_variants = 0
    nonneutral_variants = 0
    multiallelic_gt_variants = 0

    with open(vcf_path) as f:

        for line in f:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                fields = line.rstrip(
                    "\n"
                ).split("\t")

                samples = fields[9:]

                if len(samples) != len(meta):
                    raise RuntimeError(
                        f"{vcf_path}: VCF has {len(samples)} samples "
                        f"but manifest has {len(meta)} individuals"
                    )

                sample_regions = []
                sample_sites = []

                seen_indices = set()

                for sample in samples:

                    idx = sample_name_to_index(
                        sample
                    )

                    if idx in seen_indices:
                        raise RuntimeError(
                            f"{vcf_path}: duplicate VCF individual {idx}"
                        )

                    seen_indices.add(idx)

                    if idx not in meta:

                        raise RuntimeError(
                            f"{vcf_path}: VCF individual {idx} "
                            "not found in manifest"
                        )

                    region = meta[idx][
                        "region"
                    ]

                    site = meta[idx][
                        "site"
                    ]

                    sample_regions.append(
                        region
                    )

                    sample_sites.append(
                        site
                    )

                    n_individuals_region[
                        region
                    ] += 1

                continue

            if line.startswith("#"):
                continue

            if sample_regions is None:
                raise RuntimeError(
                    f"No #CHROM header found before variants in {vcf_path}"
                )

            total_variants += 1

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) < 10:
                continue

            info = parse_info(
                fields[7]
            )

            mt = info.get(
                "MT"
            )

            if mt is None:

                raise RuntimeError(
                    f"{vcf_path}: variant without MT field"
                )

            # Keep only variants for which every listed mutation
            # type is m1 (neutral).
            mt_values = mt.split(",")

            if not all(
                x == "1"
                for x in mt_values
            ):

                nonneutral_variants += 1
                continue

            neutral_variants += 1

            format_fields = fields[
                8
            ].split(":")

            try:

                gt_index = (
                    format_fields.index(
                        "GT"
                    )
                )

            except ValueError:

                raise RuntimeError(
                    f"{vcf_path}: FORMAT has no GT"
                )

            genotypes = fields[9:]

            allele_count = defaultdict(int)
            allele_number = defaultdict(int)

            weird_multiallelic = False

            for region, sample_field in zip(
                sample_regions,
                genotypes
            ):

                sample_parts = (
                    sample_field.split(":")
                )

                gt = sample_parts[
                    gt_index
                ]

                alleles = parse_gt(
                    gt
                )

                for a in alleles:

                    # Output should be biallelic because
                    # outputMultiallelics=F.
                    if a not in (
                        0,
                        1
                    ):

                        weird_multiallelic = True
                        break

                    allele_number[
                        region
                    ] += 1

                    allele_count[
                        region
                    ] += a

                if weird_multiallelic:
                    break

            if weird_multiallelic:

                multiallelic_gt_variants += 1
                neutral_variants -= 1
                continue


            # ------------------------------------------------
            # Nucleotide diversity
            #
            # contribution at one biallelic site:
            #
            # 2 * AC * (AN - AC) /
            # [ AN * (AN - 1) ]
            #
            # Sum across sites / chromosome length.
            # ------------------------------------------------

            for region in REGION_ORDER:

                an = allele_number[
                    region
                ]

                ac = allele_count[
                    region
                ]

                if an <= 1:
                    continue

                contribution = (
                    2.0 *
                    ac *
                    (an - ac) /
                    (
                        an *
                        (an - 1)
                    )
                )

                pi_sum[
                    region
                ] += contribution


            # ------------------------------------------------
            # Hudson pairwise FST
            #
            # numerator:
            # (p1-p2)^2
            # - p1(1-p1)/(n1-1)
            # - p2(1-p2)/(n2-1)
            #
            # denominator:
            # p1(1-p2) + p2(1-p1)
            #
            # Weighted FST =
            # sum(numerator) / sum(denominator)
            # ------------------------------------------------

            for i in range(
                len(REGION_ORDER)
            ):

                for j in range(
                    i + 1,
                    len(REGION_ORDER)
                ):

                    r1 = REGION_ORDER[i]
                    r2 = REGION_ORDER[j]

                    n1 = allele_number[
                        r1
                    ]

                    n2 = allele_number[
                        r2
                    ]

                    if n1 <= 1 or n2 <= 1:
                        continue

                    p1 = (
                        allele_count[r1] /
                        n1
                    )

                    p2 = (
                        allele_count[r2] /
                        n2
                    )

                    numerator = (
                        (p1 - p2) ** 2
                        -
                        (
                            p1 *
                            (1.0 - p1) /
                            (n1 - 1)
                        )
                        -
                        (
                            p2 *
                            (1.0 - p2) /
                            (n2 - 1)
                        )
                    )

                    denominator = (
                        p1 *
                        (1.0 - p2)
                        +
                        p2 *
                        (1.0 - p1)
                    )

                    if denominator <= 0:
                        continue

                    pair = (
                        r1,
                        r2
                    )

                    fst_num_sum[
                        pair
                    ] += numerator

                    fst_den_sum[
                        pair
                    ] += denominator

                    fst_sites[
                        pair
                    ] += 1


    # ========================================================
    # Validate total sample assignment
    #
    # Regional sample sizes are intentionally allowed to vary.
    # A region may contain fewer than the target sample size,
    # or may be absent entirely in a future snapshot.
    # ========================================================

    if sum(
        n_individuals_region[
            region
        ]
        for region in REGION_ORDER
    ) != len(meta):

        raise RuntimeError(
            f"{vcf_path}: regional sample counts do not sum "
            f"to manifest size {len(meta)}"
        )


    # ========================================================
    # Final pi
    # ========================================================

    pi_rows = []

    for region in REGION_ORDER:

        n_region = n_individuals_region[
            region
        ]

        if n_region == 0:

            # Population absent at the fixed regional sites.
            # Genetic diversity is not estimable, not zero.
            pi = float("nan")
            sample_status = "absent"

        else:

            pi = (
                pi_sum[
                    region
                ] /
                CHROM_LEN
            )

            if n_region < 10:
                sample_status = "low_n_lt10"
            else:
                sample_status = "ok"

        pi_rows.append({
            "scenario": scenario,
            "year": year,
            "rep": rep,
            "region": region,
            "n_individuals": n_region,
            "sample_status": sample_status,
            "chromosome_length": CHROM_LEN,
            "neutral_variants_vcf": neutral_variants,
            "pi": pi
        })


    # ========================================================
    # Final FST
    # ========================================================

    fst_rows = []

    for i in range(
        len(REGION_ORDER)
    ):

        for j in range(
            i + 1,
            len(REGION_ORDER)
        ):

            r1 = REGION_ORDER[i]
            r2 = REGION_ORDER[j]

            pair = (
                r1,
                r2
            )

            den = fst_den_sum[
                pair
            ]

            if den > 0:

                fst = (
                    fst_num_sum[
                        pair
                    ] /
                    den
                )

            else:

                fst = float(
                    "nan"
                )

            n_region1 = n_individuals_region[r1]
            n_region2 = n_individuals_region[r2]

            if n_region1 == 0 or n_region2 == 0:
                sample_status = "absent_region"
            elif n_region1 < 10 or n_region2 < 10:
                sample_status = "low_n_lt10"
            else:
                sample_status = "ok"

            fst_rows.append({
                "scenario": scenario,
                "year": year,
                "rep": rep,
                "region1": r1,
                "region2": r2,
                "pair": f"{r1}-{r2}",
                "n_individuals_region1": n_region1,
                "n_individuals_region2": n_region2,
                "sample_status": sample_status,
                "neutral_variants_vcf": neutral_variants,
                "informative_sites": fst_sites[
                    pair
                ],
                "fst_hudson": fst
            })


    qc = {
        "scenario": scenario,
        "year": year,
        "rep": rep,
        "vcf": str(vcf_path),
        "n_samples": len(meta),
        "total_variants": total_variants,
        "neutral_variants": neutral_variants,
        "nonneutral_variants": nonneutral_variants,
        "unexpected_multiallelic_gt_variants":
            multiallelic_gt_variants
    }

    return (
        pi_rows,
        fst_rows,
        qc
    )


# ============================================================
# Locate the six states used for each selected replicate
# ============================================================

vcfs = []

for rep in REPLICATES:

    required = [
        VCF_DIR / "historical" / f"{rep}_year1900.vcf",
        VCF_DIR / "historical" / f"{rep}_year2020.vcf",
        VCF_DIR / "status_quo" / f"{rep}_year2140.vcf",
        VCF_DIR / "restore_2km" / f"{rep}_year2140.vcf",
        VCF_DIR / "restore_4km" / f"{rep}_year2140.vcf",
        VCF_DIR / "restore_6km" / f"{rep}_year2140.vcf",
    ]

    for path in required:

        if not path.is_file():
            raise FileNotFoundError(
                f"Required VCF missing: {path}"
            )

        vcfs.append(path)

print(
    "Replicates:",
    len(REPLICATES)
)

print(
    "VCFs found:",
    len(vcfs)
)


# ============================================================
# Analyse
# ============================================================

all_pi = []
all_fst = []
all_qc = []

for i, vcf in enumerate(
    vcfs,
    start=1
):

    print(
        f"[{i:02d}/{len(vcfs):02d}]",
        vcf.relative_to(BASE)
    )

    pi_rows, fst_rows, qc = (
        analyse_vcf(
            vcf
        )
    )

    all_pi.extend(
        pi_rows
    )

    all_fst.extend(
        fst_rows
    )

    all_qc.append(
        qc
    )


# ============================================================
# Generic TSV writer
# ============================================================

def write_tsv(
    rows,
    path,
    fields=None
):

    if not rows:
        return

    if fields is None:
        fields = list(
            rows[0].keys()
        )

    with open(
        path,
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(
            rows
        )


# ============================================================
# Raw replicate-level results
# ============================================================

write_tsv(
    all_pi,
    OUT_DIR /
    "regional_pi_by_replicate.tsv"
)

write_tsv(
    all_fst,
    OUT_DIR /
    "regional_fst_by_replicate.tsv"
)

write_tsv(
    all_qc,
    OUT_DIR /
    "regional_vcf_qc.tsv"
)


# ============================================================
# Summaries across five replicates
# ============================================================

pi_groups = defaultdict(
    list
)

for row in all_pi:

    key = (
        row["scenario"],
        row["year"],
        row["region"]
    )

    if not math.isnan(
        row["pi"]
    ):

        pi_groups[
            key
        ].append(
            row["pi"]
        )


pi_summary = []

for key, values in sorted(
    pi_groups.items()
):

    scenario, year, region = key

    pi_summary.append({
        "scenario": scenario,
        "year": year,
        "region": region,
        "n_reps": len(values),
        "mean_pi": mean(values),
        "sd_pi": sample_sd(values),
        "min_pi": min(values),
        "max_pi": max(values)
    })


write_tsv(
    pi_summary,
    OUT_DIR /
    "regional_pi_summary.tsv"
)


fst_groups = defaultdict(
    list
)

for row in all_fst:

    key = (
        row["scenario"],
        row["year"],
        row["region1"],
        row["region2"],
        row["pair"]
    )

    if not math.isnan(
        row["fst_hudson"]
    ):

        fst_groups[
            key
        ].append(
            row["fst_hudson"]
        )


fst_summary = []

for key, values in sorted(
    fst_groups.items()
):

    (
        scenario,
        year,
        r1,
        r2,
        pair
    ) = key

    fst_summary.append({
        "scenario": scenario,
        "year": year,
        "region1": r1,
        "region2": r2,
        "pair": pair,
        "n_reps": len(values),
        "mean_fst": mean(values),
        "sd_fst": sample_sd(values),
        "min_fst": min(values),
        "max_fst": max(values)
    })


write_tsv(
    fst_summary,
    OUT_DIR /
    "regional_fst_summary.tsv"
)


# ============================================================
# Dictionaries for paired comparisons
# ============================================================

pi_lookup = {
    (
        r["scenario"],
        r["year"],
        r["rep"],
        r["region"]
    ):
    r["pi"]

    for r in all_pi
}


fst_lookup = {
    (
        r["scenario"],
        r["year"],
        r["rep"],
        r["pair"]
    ):
    r["fst_hudson"]

    for r in all_fst
}


# ============================================================
# PI: paired temporal changes
# ============================================================

pi_changes = []

REPS = [
    f"rep{i:03d}"
    for i in range(
        1,
        6
    )
]

FUTURE_SCENARIOS = [
    "status_quo",
    "restore_2km",
    "restore_4km",
    "restore_6km"
]


for rep in REPS:

    for region in REGION_ORDER:

        p1900 = pi_lookup[
            (
                "shared_history",
                1900,
                rep,
                region
            )
        ]

        p2020 = pi_lookup[
            (
                "shared_history",
                2020,
                rep,
                region
            )
        ]

        pi_changes.append({
            "comparison":
                "1900_to_2020",
            "scenario":
                "shared_history",
            "rep":
                rep,
            "region":
                region,
            "start_pi":
                p1900,
            "end_pi":
                p2020,
            "delta_pi":
                p2020 - p1900,
            "retention":
                safe_ratio(
                    p2020,
                    p1900
                ),
            "percent_change":
                100.0 *
                (
                    safe_ratio(
                        p2020,
                        p1900
                    ) - 1.0
                )
        })


        for scenario in FUTURE_SCENARIOS:

            p2140 = pi_lookup[
                (
                    scenario,
                    2140,
                    rep,
                    region
                )
            ]

            pi_changes.append({
                "comparison":
                    "2020_to_2140",
                "scenario":
                    scenario,
                "rep":
                    rep,
                "region":
                    region,
                "start_pi":
                    p2020,
                "end_pi":
                    p2140,
                "delta_pi":
                    p2140 - p2020,
                "retention":
                    safe_ratio(
                        p2140,
                        p2020
                    ),
                "percent_change":
                    100.0 *
                    (
                        safe_ratio(
                            p2140,
                            p2020
                        ) - 1.0
                    )
            })


        sq2140 = pi_lookup[
            (
                "status_quo",
                2140,
                rep,
                region
            )
        ]

        for scenario in [
            "restore_2km",
            "restore_4km",
            "restore_6km"
        ]:

            rest2140 = pi_lookup[
                (
                    scenario,
                    2140,
                    rep,
                    region
                )
            ]

            pi_changes.append({
                "comparison":
                    "2140_restoration_vs_statusquo",
                "scenario":
                    scenario,
                "rep":
                    rep,
                "region":
                    region,
                "start_pi":
                    sq2140,
                "end_pi":
                    rest2140,
                "delta_pi":
                    rest2140 - sq2140,
                "retention":
                    safe_ratio(
                        rest2140,
                        sq2140
                    ),
                "percent_change":
                    100.0 *
                    (
                        safe_ratio(
                            rest2140,
                            sq2140
                        ) - 1.0
                    )
            })


write_tsv(
    pi_changes,
    OUT_DIR /
    "regional_pi_paired_changes.tsv"
)


# ============================================================
# FST: paired temporal changes
# ============================================================

PAIRS = [
    f"{REGION_ORDER[i]}-{REGION_ORDER[j]}"
    for i in range(
        len(REGION_ORDER)
    )
    for j in range(
        i + 1,
        len(REGION_ORDER)
    )
]


fst_changes = []

for rep in REPS:

    for pair in PAIRS:

        f1900 = fst_lookup[
            (
                "shared_history",
                1900,
                rep,
                pair
            )
        ]

        f2020 = fst_lookup[
            (
                "shared_history",
                2020,
                rep,
                pair
            )
        ]

        fst_changes.append({
            "comparison":
                "1900_to_2020",
            "scenario":
                "shared_history",
            "rep":
                rep,
            "pair":
                pair,
            "start_fst":
                f1900,
            "end_fst":
                f2020,
            "delta_fst":
                f2020 - f1900
        })


        for scenario in FUTURE_SCENARIOS:

            f2140 = fst_lookup[
                (
                    scenario,
                    2140,
                    rep,
                    pair
                )
            ]

            fst_changes.append({
                "comparison":
                    "2020_to_2140",
                "scenario":
                    scenario,
                "rep":
                    rep,
                "pair":
                    pair,
                "start_fst":
                    f2020,
                "end_fst":
                    f2140,
                "delta_fst":
                    f2140 - f2020
            })


        sq2140 = fst_lookup[
            (
                "status_quo",
                2140,
                rep,
                pair
            )
        ]

        for scenario in [
            "restore_2km",
            "restore_4km",
            "restore_6km"
        ]:

            rest2140 = fst_lookup[
                (
                    scenario,
                    2140,
                    rep,
                    pair
                )
            ]

            fst_changes.append({
                "comparison":
                    "2140_restoration_vs_statusquo",
                "scenario":
                    scenario,
                "rep":
                    rep,
                "pair":
                    pair,
                "start_fst":
                    sq2140,
                "end_fst":
                    rest2140,
                "delta_fst":
                    rest2140 - sq2140
            })


write_tsv(
    fst_changes,
    OUT_DIR /
    "regional_fst_paired_changes.tsv"
)


# ============================================================
# Summarise paired changes
# ============================================================

def summarise_changes(
    rows,
    group_fields,
    value_field
):

    groups = defaultdict(
        list
    )

    for row in rows:

        key = tuple(
            row[x]
            for x in group_fields
        )

        value = row[
            value_field
        ]

        if not math.isnan(
            value
        ):

            groups[
                key
            ].append(
                value
            )

    result = []

    for key, values in sorted(
        groups.items()
    ):

        row = dict(
            zip(
                group_fields,
                key
            )
        )

        row.update({
            "n_reps":
                len(values),
            "mean":
                mean(values),
            "sd":
                sample_sd(values),
            "min":
                min(values),
            "max":
                max(values)
        })

        result.append(
            row
        )

    return result


pi_change_summary = (
    summarise_changes(
        pi_changes,
        [
            "comparison",
            "scenario",
            "region"
        ],
        "percent_change"
    )
)


write_tsv(
    pi_change_summary,
    OUT_DIR /
    "regional_pi_change_summary.tsv"
)


fst_change_summary = (
    summarise_changes(
        fst_changes,
        [
            "comparison",
            "scenario",
            "pair"
        ],
        "delta_fst"
    )
)


write_tsv(
    fst_change_summary,
    OUT_DIR /
    "regional_fst_change_summary.tsv"
)


# ============================================================
# Terminal summaries
# ============================================================

print()
print(
    "============================================"
)
print(
    "VCF QC"
)
print(
    "============================================"
)

for row in all_qc:

    print(
        row["scenario"],
        row["year"],
        row["rep"],
        "total=",
        row["total_variants"],
        "neutral=",
        row["neutral_variants"],
        "selected=",
        row["nonneutral_variants"],
        "weird=",
        row[
            "unexpected_multiallelic_gt_variants"
        ]
    )


print()
print(
    "============================================"
)
print(
    "MEAN REGIONAL PI"
)
print(
    "============================================"
)

for row in pi_summary:

    print(
        f"{row['scenario']:15s}",
        row["year"],
        f"{row['region']:>2s}",
        f"pi={row['mean_pi']:.8f}",
        f"sd={row['sd_pi']:.8f}"
    )


print()
print(
    "============================================"
)
print(
    "MEAN REGIONAL HUDSON FST"
)
print(
    "============================================"
)

for row in fst_summary:

    print(
        f"{row['scenario']:15s}",
        row["year"],
        f"{row['pair']:>5s}",
        f"FST={row['mean_fst']:.6f}",
        f"sd={row['sd_fst']:.6f}"
    )


print()
print(
    "============================================"
)
print(
    "DONE"
)
print(
    "============================================"
)

print(
    "Output directory:",
    OUT_DIR
)
