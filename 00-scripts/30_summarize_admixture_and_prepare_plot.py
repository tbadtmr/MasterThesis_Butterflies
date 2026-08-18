#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import re

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "tabea_work" / "08_population_analysis_final" / "skane_only"

ADMIX = BASE / "combined" / "results" / "structure" / "admixture"
PCA   = BASE / "combined" / "results" / "structure" / "pca"

OUT = BASE / "combined" / "results" / "structure" / "final_plot"
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. PARSE ALL ADMIXTURE RUNS
# ============================================================

records = []

for repdir in sorted(ADMIX.glob("rep*")):

    rep = repdir.name

    for kdir in sorted(
        repdir.glob("K*"),
        key=lambda p: int(p.name[1:])
    ):

        K = int(kdir.name[1:])

        for rundir in sorted(kdir.glob("run*")):

            logfile = rundir / "run.log"

            if not logfile.exists():
                continue

            txt = logfile.read_text(
                errors="ignore"
            )

            cv_match = re.search(
                r"CV error\s*\(K\s*=\s*\d+\)\s*:\s*"
                r"([-+0-9.eE]+)",
                txt
            )

            ll_matches = re.findall(
                r"Loglikelihood\s*:\s*([-+0-9.eE]+)",
                txt
            )

            cv = (
                float(cv_match.group(1))
                if cv_match
                else np.nan
            )

            ll = (
                float(ll_matches[-1])
                if ll_matches
                else np.nan
            )

            qfiles = list(rundir.glob("*.Q"))

            qfile = (
                str(qfiles[0])
                if qfiles
                else None
            )

            records.append({
                "replicate": rep,
                "K": K,
                "run": rundir.name,
                "cv_error": cv,
                "loglikelihood": ll,
                "q_file": qfile,
                "run_dir": str(rundir)
            })

runs = pd.DataFrame(records)

if runs.empty:
    raise SystemExit(
        "ERROR: no ADMIXTURE runs found"
    )

runs.to_csv(
    OUT / "admixture_all_runs.tsv",
    sep="\t",
    index=False
)

print("ADMIXTURE runs parsed:", len(runs))

# ============================================================
# 2. SUMMARY WITHIN EACH REPLICATE x K
# ============================================================

rep_rows = []

for (rep, K), g in runs.groupby(
    ["replicate", "K"]
):

    good_ll = g[
        np.isfinite(g["loglikelihood"])
    ].sort_values(
        "loglikelihood",
        ascending=False
    )

    top3 = good_ll.head(3)

    if len(top3) >= 3:
        spread = (
            top3["loglikelihood"].max()
            - top3["loglikelihood"].min()
        )
    else:
        spread = np.nan

    best = (
        good_ll.iloc[0]
        if len(good_ll)
        else g.iloc[0]
    )

    rep_rows.append({
        "replicate": rep,
        "K": K,

        "n_runs": len(g),

        "mean_cv": g["cv_error"].mean(),
        "sd_cv": g["cv_error"].std(),

        "best_loglikelihood":
            best["loglikelihood"],

        "top3_ll_spread":
            spread,

        # Zach-style convergence criterion:
        # top three likelihoods within 2 LL units
        "top3_converged":
            bool(
                np.isfinite(spread)
                and spread <= 2
            ),

        "best_run":
            best["run"],

        "best_q_file":
            best["q_file"],

        "best_run_dir":
            best["run_dir"]
    })

rep_summary = pd.DataFrame(rep_rows)

rep_summary.to_csv(
    OUT / "admixture_by_replicate_K.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 3. SUMMARY ACROSS SIMULATION REPLICATES
# ============================================================

k_summary = (
    rep_summary
    .groupby("K")
    .agg(
        n_replicates=("replicate", "nunique"),
        mean_cv=("mean_cv", "mean"),
        sd_cv=("mean_cv", "std"),
        min_cv=("mean_cv", "min"),
        max_cv=("mean_cv", "max"),
        converged_replicates=(
            "top3_converged",
            "sum"
        ),
        mean_top3_spread=(
            "top3_ll_spread",
            "mean"
        )
    )
    .reset_index()
)

# Primary ADMIXTURE K choice:
# lowest mean CV across independent simulation replicates.
# Lowest mean cross-validation error.
CV_BEST_K = int(
    k_summary
    .sort_values(["mean_cv", "K"])
    .iloc[0]["K"]
)

# For visualization, retain the highest K for which all
# independent simulation replicates satisfy the top-three
# likelihood convergence criterion.
fully_converged = k_summary[
    k_summary["converged_replicates"] ==
    k_summary["n_replicates"]
]

BEST_K = int(
    fully_converged["K"].max()
)

k_summary["best_by_cv"] = (
    k_summary["K"] == CV_BEST_K
)

k_summary["selected_for_plot"] = (
    k_summary["K"] == BEST_K
)

k_summary.to_csv(
    OUT / "admixture_K_summary.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 4. REPRESENTATIVE SIMULATION REPLICATE
#
# Prefer converged runs, then choose replicate whose CV is
# closest to the median CV at the selected K.
# ============================================================

bestk = rep_summary[
    rep_summary["K"] == BEST_K
].copy()

conv = bestk[
    bestk["top3_converged"]
].copy()

candidate = (
    conv
    if len(conv) > 0
    else bestk
)

median_cv = candidate["mean_cv"].median()

candidate["distance_to_median_cv"] = (
    candidate["mean_cv"] - median_cv
).abs()

chosen = candidate.sort_values(
    [
        "distance_to_median_cv",
        "top3_ll_spread",
        "replicate"
    ]
).iloc[0]

REP = chosen["replicate"]
BEST_RUN = chosen["best_run"]
QFILE = Path(chosen["best_q_file"])
RUNDIR = Path(chosen["best_run_dir"])

print()
print("============================================")
print("ADMIXTURE MODEL SELECTION")
print("============================================")
print()
print(k_summary.to_string(index=False))
print()
print("Lowest-CV K       :", CV_BEST_K)
print("Selected K         :", BEST_K)
print("Representative rep:", REP)
print("Selected run      :", BEST_RUN)
print("Q file            :", QFILE)
print()

# ============================================================
# 5. READ THE SIX MANIFESTS FOR THIS REPLICATE
# ============================================================

states = [
    ("historical",   1900, "1900"),
    ("historical",   2020, "2020"),
    ("status_quo",   2140, "SQ"),
    ("restore_2km",  2140, "R2"),
    ("restore_4km",  2140, "R4"),
    ("restore_6km",  2140, "R6"),
]

meta_frames = []

for scenario, year, state in states:

    f = (
        BASE /
        "combined" /
        "sampling" /
        "manifests" /
        scenario /
        f"{REP}_year{year}_manifest.tsv"
    )

    if not f.exists():
        raise SystemExit(
            f"ERROR: manifest missing: {f}"
        )

    m = pd.read_csv(
        f,
        sep="\t"
    )

    m["state"] = state
    m["scenario_plot"] = scenario

    meta_frames.append(m)

meta = pd.concat(
    meta_frames,
    ignore_index=True
)

# Make sure replicate naming is standardized
if "replicate" not in meta.columns:
    if "rep" in meta.columns:
        meta["replicate"] = meta["rep"]

# ============================================================
# 6. FIND ADMIXTURE .FAM FILE
# ============================================================

fam_candidates = list(
    RUNDIR.glob("*.fam")
)

if not fam_candidates:

    fam_candidates = list(
        (ADMIX / REP).rglob("*.fam")
    )

if not fam_candidates:

    fam_candidates = list(
        PCA.glob(f"{REP}*.fam")
    )

if not fam_candidates:

    fam_candidates = list(
        PCA.rglob(f"*{REP}*.fam")
    )

if not fam_candidates:
    raise SystemExit(
        "ERROR: could not locate the PLINK .fam used by ADMIXTURE"
    )

FAM = fam_candidates[0]

fam = pd.read_csv(
    FAM,
    sep=r"\s+",
    header=None
)

fam.columns = [
    "FID",
    "IID",
    "father",
    "mother",
    "sex_plink",
    "phenotype"
]

print("FAM:", FAM)
print("FAM individuals:", len(fam))

# ============================================================
# 7. MATCH FAM/PCA IDS TO MANIFEST
# ============================================================

def infer_record(iid):

    s = str(iid)

    # ---------------------------------------------
    # A. Exact manifest sample_id
    # ---------------------------------------------

    if "sample_id" in meta.columns:

        hit = meta[
            meta["sample_id"].astype(str) == s
        ]

        if len(hit) == 1:
            return hit.iloc[0]

    # ---------------------------------------------
    # B. Manifest sample_id contained in joint ID
    # ---------------------------------------------

    if "sample_id" in meta.columns:

        hits = meta[
            meta["sample_id"]
            .astype(str)
            .apply(lambda x: x in s)
        ]

        if len(hits) == 1:
            return hits.iloc[0]

    # ---------------------------------------------
    # C. Match simulation ind_index
    # ---------------------------------------------

    mm = re.search(
        r"i(\d+)",
        s
    )

    if mm:

        idx = int(mm.group(1))

        hits = meta[
            meta["ind_index"] == idx
        ].copy()

        if len(hits) == 1:
            return hits.iloc[0]

        if len(hits) > 1:

            # Use state/year/scenario hints in joint ID
            sl = s.lower()

            state_hint = None

            if "1900" in sl:
                state_hint = "1900"

            elif "2020" in sl:
                state_hint = "2020"

            elif (
                "status" in sl or
                "status_quo" in sl or
                "sq2140" in sl or
                "_sq" in sl
            ):
                state_hint = "SQ"

            elif (
                "restore_2" in sl or
                "_r2" in sl or
                sl.startswith("r2_")
            ):
                state_hint = "R2"

            elif (
                "restore_4" in sl or
                "_r4" in sl or
                sl.startswith("r4_")
            ):
                state_hint = "R4"

            elif (
                "restore_6" in sl or
                "_r6" in sl or
                sl.startswith("r6_")
            ):
                state_hint = "R6"

            if state_hint is not None:

                hh = hits[
                    hits["state"] == state_hint
                ]

                if len(hh) == 1:
                    return hh.iloc[0]

    return None


matched = []

unmatched = []

for _, frow in fam.iterrows():

    rec = infer_record(
        frow["IID"]
    )

    if rec is None:

        unmatched.append(
            str(frow["IID"])
        )

        matched.append({
            "FID": frow["FID"],
            "IID": frow["IID"]
        })

    else:

        d = rec.to_dict()

        d["FID"] = frow["FID"]
        d["IID"] = frow["IID"]

        matched.append(d)

plot_meta = pd.DataFrame(
    matched
)

if unmatched:

    pd.DataFrame({
        "IID": unmatched
    }).to_csv(
        OUT / "UNMATCHED_SAMPLE_IDS.tsv",
        sep="\t",
        index=False
    )

    raise SystemExit(
        f"ERROR: {len(unmatched)} FAM IDs could not be matched. "
        f"See {OUT/'UNMATCHED_SAMPLE_IDS.tsv'}"
    )

plot_meta["fam_order"] = np.arange(
    1,
    len(plot_meta) + 1
)

# ============================================================
# 8. ATTACH SELECTED Q MATRIX
# ============================================================

Q = pd.read_csv(
    QFILE,
    sep=r"\s+",
    header=None
)

if len(Q) != len(plot_meta):

    raise SystemExit(
        f"ERROR: Q rows ({len(Q)}) != FAM rows ({len(plot_meta)})"
    )

Q.columns = [
    f"Q{i}"
    for i in range(
        1,
        BEST_K + 1
    )
]

qdata = pd.concat(
    [
        plot_meta.reset_index(drop=True),
        Q.reset_index(drop=True)
    ],
    axis=1
)

qdata.to_csv(
    OUT / "representative_admixture_Q.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 9. PCA DATA FOR SAME REPRESENTATIVE REPLICATE
# ============================================================

eigvec = (
    PCA /
    f"{REP}_six_states_PCA.eigenvec"
)

eigval = (
    PCA /
    f"{REP}_six_states_PCA.eigenval"
)

if not eigvec.exists():

    raise SystemExit(
        f"ERROR: PCA eigenvec missing: {eigvec}"
    )

pc = pd.read_csv(
    eigvec,
    sep=r"\s+",
    header=None
)

pc.columns = (
    ["FID", "IID"] +
    [
        f"PC{i}"
        for i in range(
            1,
            pc.shape[1] - 1
        )
    ]
)

# Match PCA IDs using metadata prepared from fam.
lookup = qdata[
    [
        "IID",
        "state",
        "region",
        "site",
        "site_number",
        "ind_index"
    ]
].drop_duplicates()

pc = pc.merge(
    lookup,
    on="IID",
    how="left"
)

if pc["state"].isna().any():

    bad = pc.loc[
        pc["state"].isna(),
        "IID"
    ]

    raise SystemExit(
        "ERROR: PCA IDs failed metadata match: "
        + ", ".join(
            bad.astype(str).head(10)
        )
    )

pc.to_csv(
    OUT / "representative_PCA.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 10. SELECTION RECORD
# ============================================================

with open(
    OUT / "admixture_selection.txt",
    "w"
) as fh:

    fh.write(
        f"best_K_by_mean_CV\t{BEST_K}\n"
    )

    fh.write(
        f"representative_replicate\t{REP}\n"
    )

    fh.write(
        f"selected_run\t{BEST_RUN}\n"
    )

    fh.write(
        f"mean_CV_representative\t{chosen['mean_cv']}\n"
    )

    fh.write(
        f"top3_LL_spread\t{chosen['top3_ll_spread']}\n"
    )

    fh.write(
        f"top3_converged\t{chosen['top3_converged']}\n"
    )

print("Prepared plotting data in:")
print(OUT)
print()
print("DONE")
