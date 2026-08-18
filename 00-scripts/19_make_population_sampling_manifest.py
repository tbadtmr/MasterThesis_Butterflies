#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict

# ============================================================
# Fixed regional sampling
#
# Regions and sites are read from the supplied design file.
# The design must contain 3 fixed sites per region.
#
# radius = 4 model units
# target = up to 20 individuals / site
#
# => up to 60 individuals / region
#
# Sampling is deterministic using SHA256 ranking.
# ============================================================

if len(sys.argv) != 4:
    sys.exit(
        "Usage: python3 19_make_population_sampling_manifest.py "
        "WORK_DIR DESIGN_FILE CLUSTER"
    )

WORK = Path(sys.argv[1]).expanduser().resolve()
DESIGN = Path(sys.argv[2]).expanduser().resolve()
CLUSTER = sys.argv[3]

POSITIONS = WORK / "positions"
SAMPLING = WORK / "sampling"
MANIFEST_DIR = SAMPLING / "manifests"

RADIUS = 4.0
TARGET_N = 20
MASTER_SEED = os.environ.get(
    "POP_SAMPLING_SEED",
    "20260817_skane_structure"
)

SAMPLING.mkdir(parents=True, exist_ok=True)
MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

REGION_ORDER = []
SCENARIO_ORDER = [
    "historical",
    "status_quo",
    "restore_2km",
    "restore_4km",
    "restore_6km",
]

# ------------------------------------------------------------
# Read fixed regional sampling design
# ------------------------------------------------------------

sites = []
region_counter = defaultdict(int)

with DESIGN.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    required = {"site", "region", "x0", "y0"}

    if not required.issubset(reader.fieldnames or []):
        raise RuntimeError(
            f"Design file lacks required columns: {required}"
        )

    for row in reader:
        region = row["region"]

        if region not in REGION_ORDER:
            REGION_ORDER.append(region)

        region_counter[region] += 1

        sites.append({
            "site": row["site"],
            "region": region,
            "site_number": region_counter[region],
            "x0": float(row["x0"]),
            "y0": float(row["y0"]),
        })

if not sites:
    raise RuntimeError("Sampling design contains no sites")

for region in REGION_ORDER:
    if region_counter[region] != 3:
        raise RuntimeError(
            f"{region}: expected 3 sites; "
            f"found {region_counter[region]}"
        )

print("Sampling design")
print("----------------")
for s in sites:
    print(
        s["region"],
        s["site_number"],
        s["site"],
        s["x0"],
        s["y0"]
    )

print()
print("Radius:", RADIUS)
print("Target/site:", TARGET_N)
print()

# ------------------------------------------------------------
# Deterministic pseudo-random rank
# ------------------------------------------------------------

def rank_individual(scenario, rep, year, site, ind_index):

    key = (
        f"{MASTER_SEED}|{scenario}|{rep}|"
        f"{year}|{site}|{ind_index}"
    )

    return hashlib.sha256(
        key.encode("utf-8")
    ).hexdigest()

# ------------------------------------------------------------
# Find position files
# ------------------------------------------------------------

position_files = []

for scenario in SCENARIO_ORDER:

    scenario_dir = POSITIONS / scenario

    if not scenario_dir.exists():
        continue

    for f in scenario_dir.glob("*_positions.tsv"):

        if f.name.startswith("test_"):
            continue

        m = re.fullmatch(
            r"(rep\d+)_year(\d+)_positions\.tsv",
            f.name
        )

        if m is None:
            continue

        position_files.append(
            (
                scenario,
                m.group(1),
                int(m.group(2)),
                f
            )
        )

position_files.sort(
    key=lambda x: (
        SCENARIO_ORDER.index(x[0]),
        x[1],
        x[2]
    )
)

if not position_files:
    raise RuntimeError(
        f"No position files found under {POSITIONS}"
    )

print("Position files:", len(position_files))
print()

manifest_fields = [
    "sample_id",
    "cluster",
    "scenario",
    "replicate",
    "year",
    "region",
    "site",
    "site_number",
    "site_x",
    "site_y",
    "radius_model_units",
    "target_n",
    "n_available_at_site",
    "ind_index",
    "sex",
    "age",
    "x",
    "y",
    "distance_to_site",
]

availability_fields = [
    "cluster",
    "scenario",
    "replicate",
    "year",
    "region",
    "site",
    "site_number",
    "site_x",
    "site_y",
    "radius_model_units",
    "target_n",
    "n_available",
    "n_selected",
    "target_met",
]

all_selected = []
all_availability = []
snapshot_summary = []
region_summary = []

# ------------------------------------------------------------
# Process snapshots
# ------------------------------------------------------------

for scenario, rep, year, pos_file in position_files:

    print(
        f"Sampling {scenario:15s} "
        f"{rep} {year}"
    )

    individuals = []

    with pos_file.open() as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        required = {
            "ind_index",
            "sex",
            "age",
            "x",
            "y",
        }

        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError(
                f"{pos_file}: missing required columns"
            )

        for row in reader:

            individuals.append({
                "ind_index": row["ind_index"],
                "sex": row["sex"],
                "age": row["age"],
                "x": float(row["x"]),
                "y": float(row["y"]),
            })

    selected_snapshot = []
    availability_snapshot = []

    for s in sites:

        candidates = []

        for ind in individuals:

            dx = ind["x"] - s["x0"]
            dy = ind["y"] - s["y0"]

            dist2 = dx * dx + dy * dy

            if dist2 <= RADIUS * RADIUS:

                distance = dist2 ** 0.5

                rank = rank_individual(
                    scenario,
                    rep,
                    year,
                    s["site"],
                    ind["ind_index"],
                )

                candidates.append(
                    (rank, ind, distance)
                )

        candidates.sort(key=lambda z: z[0])

        n_available = len(candidates)
        n_selected = min(
            TARGET_N,
            n_available
        )

        availability_row = {
            "cluster": CLUSTER,
            "scenario": scenario,
            "replicate": rep,
            "year": year,
            "region": s["region"],
            "site": s["site"],
            "site_number": s["site_number"],
            "site_x": s["x0"],
            "site_y": s["y0"],
            "radius_model_units": RADIUS,
            "target_n": TARGET_N,
            "n_available": n_available,
            "n_selected": n_selected,
            "target_met": n_available >= TARGET_N,
        }

        availability_snapshot.append(
            availability_row
        )

        all_availability.append(
            availability_row
        )

        for _, ind, distance in candidates[:TARGET_N]:

            sample_id = (
                f"{scenario}_{rep}_{year}_"
                f"{s['region']}_{s['site']}_"
                f"i{ind['ind_index']}"
            )

            row = {
                "sample_id": sample_id,
                "cluster": CLUSTER,
                "scenario": scenario,
                "replicate": rep,
                "year": year,
                "region": s["region"],
                "site": s["site"],
                "site_number": s["site_number"],
                "site_x": s["x0"],
                "site_y": s["y0"],
                "radius_model_units": RADIUS,
                "target_n": TARGET_N,
                "n_available_at_site": n_available,
                "ind_index": ind["ind_index"],
                "sex": ind["sex"],
                "age": ind["age"],
                "x": ind["x"],
                "y": ind["y"],
                "distance_to_site": distance,
            }

            selected_snapshot.append(row)
            all_selected.append(row)

    # The circles should not overlap enough to sample
    # the same individual twice.
    ids = [
        r["ind_index"]
        for r in selected_snapshot
    ]

    if len(ids) != len(set(ids)):
        raise RuntimeError(
            f"Duplicate sampled individual across sites: "
            f"{scenario} {rep} {year}"
        )

    # --------------------------------------------------------
    # Per-snapshot manifest
    # --------------------------------------------------------

    out_dir = MANIFEST_DIR / scenario
    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    outfile = (
        out_dir /
        f"{rep}_year{year}_manifest.tsv"
    )

    with outfile.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=manifest_fields,
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(
            selected_snapshot
        )

    # --------------------------------------------------------
    # Snapshot summary
    # --------------------------------------------------------

    site_n = [
        r["n_selected"]
        for r in availability_snapshot
    ]

    snapshot_summary.append({
        "cluster": CLUSTER,
        "scenario": scenario,
        "replicate": rep,
        "year": year,
        "n_sites": len(sites),
        "n_sites_target_met": sum(
            1
            for r in availability_snapshot
            if r["target_met"]
        ),
        "n_individuals_selected":
            len(selected_snapshot),
        "min_selected_per_site":
            min(site_n),
        "max_selected_per_site":
            max(site_n),
    })

    # --------------------------------------------------------
    # Region summary
    # --------------------------------------------------------

    for region in REGION_ORDER:

        rows_region = [
            r
            for r in selected_snapshot
            if r["region"] == region
        ]

        region_summary.append({
            "cluster": CLUSTER,
            "scenario": scenario,
            "replicate": rep,
            "year": year,
            "region": region,
            "n_selected":
                len(rows_region),
        })

# ============================================================
# Write combined outputs
# ============================================================

def write_table(path, rows):

    if not rows:
        return

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=list(rows[0].keys()),
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)

write_table(
    SAMPLING /
    "population_sampling_manifest.tsv",
    all_selected
)

write_table(
    SAMPLING /
    "site_availability.tsv",
    all_availability
)

write_table(
    SAMPLING /
    "snapshot_sampling_summary.tsv",
    snapshot_summary
)

write_table(
    SAMPLING /
    "region_sampling_summary.tsv",
    region_summary
)

print()
print("====================================")
print("DONE")
print("Snapshots:", len(position_files))
print("Selected individuals:", len(all_selected))
print(
    "Site/snapshot combinations below 20:",
    sum(
        1
        for r in all_availability
        if not r["target_met"]
    )
)
print("====================================")
