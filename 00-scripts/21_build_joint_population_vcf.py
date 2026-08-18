#!/usr/bin/env python3

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument(
    "--output",
    required=True
)

parser.add_argument(
    "--input",
    action="append",
    required=True,
    help="LABEL=VCF"
)

args = parser.parse_args()


# ============================================================
# Parse inputs
# ============================================================

inputs = []

for item in args.input:

    if "=" not in item:
        raise SystemExit(
            f"Bad --input: {item}"
        )

    label, path = item.split(
        "=",
        1
    )

    inputs.append(
        (
            label,
            Path(path)
        )
    )


for label, path in inputs:

    if not path.exists():
        raise SystemExit(
            f"Missing VCF: {path}"
        )


# ============================================================
# Helpers
# ============================================================

def neutral_record(fields):

    if len(fields) < 10:
        return False

    if "," in fields[4]:
        return False

    mt = None

    for item in fields[7].split(";"):

        if item.startswith("MT="):
            mt = item.split(
                "=",
                1
            )[1]

    if mt is None:
        return False

    # Only m1 / MT1 neutral mutations
    return all(
        x == "1"
        for x in mt.split(",")
    )


def get_gt(sample_field, format_fields):

    parts = sample_field.split(":")

    try:
        i = format_fields.index("GT")
    except ValueError:
        raise RuntimeError(
            "FORMAT does not contain GT"
        )

    return parts[i]


def scan_keys(path):

    keys = set()

    n_neutral = 0

    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if not neutral_record(fields):
                continue

            n_neutral += 1

            key = (
                fields[0],
                int(fields[1]),
                fields[3],
                fields[4]
            )

            keys.add(
                key
            )

    return keys, n_neutral


def read_samples(path):

    with open(path) as f:

        for line in f:

            if line.startswith("#CHROM"):

                fields = line.rstrip(
                    "\n"
                ).split("\t")

                return fields[9:]

    raise RuntimeError(
        f"No #CHROM line in {path}"
    )


# ============================================================
# Find variants shared by every state
# ============================================================

common = None

for label, path in inputs:

    keys, n = scan_keys(path)

    print(
        label,
        "neutral variants =",
        n
    )

    if common is None:

        common = keys

    else:

        common &= keys

    print(
        "  current shared variants =",
        len(common)
    )


if not common:

    raise RuntimeError(
        "No neutral variants shared across inputs."
    )


# ============================================================
# Sample names
# ============================================================

all_sample_names = []

sample_lists = {}

for label, path in inputs:

    samples = read_samples(
        path
    )

    sample_lists[label] = samples

    prefixed = [
        f"{label}__{x}"
        for x in samples
    ]

    all_sample_names.extend(
        prefixed
    )


if len(all_sample_names) != len(set(all_sample_names)):

    raise RuntimeError(
        "Combined sample IDs are not unique."
    )


# ============================================================
# Read genotypes for common variants
# ============================================================

geno_by_input = {}

metadata = {}
order = None


for file_number, (label, path) in enumerate(inputs):

    dat = {}

    this_order = []

    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if not neutral_record(fields):
                continue

            key = (
                fields[0],
                int(fields[1]),
                fields[3],
                fields[4]
            )

            if key not in common:
                continue

            fmt = fields[8].split(":")

            gts = [
                get_gt(
                    x,
                    fmt
                )
                for x in fields[9:]
            ]

            if len(gts) != len(
                sample_lists[label]
            ):

                raise RuntimeError(
                    f"Sample count mismatch: {label}"
                )

            dat[key] = "\t".join(
                gts
            )

            this_order.append(
                key
            )

            if file_number == 0:

                metadata[key] = (
                    fields[0],
                    fields[1],
                    f"{fields[0]}:{fields[1]}",
                    fields[3],
                    fields[4]
                )

    geno_by_input[
        label
    ] = dat

    if file_number == 0:

        order = [
            key
            for key in this_order
            if key in common
        ]


# ============================================================
# Validation
# ============================================================

for label, dat in geno_by_input.items():

    missing = common - set(
        dat.keys()
    )

    if missing:

        raise RuntimeError(
            f"{label}: missing {len(missing)} common variants"
        )


# ============================================================
# Write merged VCF
# ============================================================

out = Path(
    args.output
)

out.parent.mkdir(
    parents=True,
    exist_ok=True
)


with open(out, "w") as f:

    f.write(
        "##fileformat=VCFv4.2\n"
    )

    f.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )

    header = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT"
    ] + all_sample_names

    f.write(
        "\t".join(header)
        + "\n"
    )

    written = 0

    for key in order:

        if key not in common:
            continue

        chrom, pos, vid, ref, alt = (
            metadata[key]
        )

        row = [
            chrom,
            pos,
            vid,
            ref,
            alt,
            ".",
            "PASS",
            ".",
            "GT"
        ]

        for label, path in inputs:

            row.append(
                geno_by_input[label][key]
            )

        f.write(
            "\t".join(row)
            + "\n"
        )

        written += 1


print()
print(
    "========================================"
)
print(
    "Joint VCF written:",
    out
)
print(
    "Input states:",
    len(inputs)
)
print(
    "Samples:",
    len(all_sample_names)
)
print(
    "Shared neutral variants:",
    written
)
print(
    "========================================"
)
