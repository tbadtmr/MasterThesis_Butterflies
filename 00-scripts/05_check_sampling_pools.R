# ============================================================
# Check locality-matched sampling pools
# Cyaniris semiargus empirical validation
#
# Uses the exact coordinates/sample sizes in validation_targets.tsv.
#
# Important:
#   - all individuals are age 0 in these snapshots, so no age filter
#   - skane_only excludes NSmaland targets
#   - radius is read from validation_targets.tsv
# ============================================================

options(stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Read and validate targets
# ------------------------------------------------------------

targets <- read.delim(
    "validation_targets.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_target_cols <- c(
    "target",
    "region",
    "empirical_year",
    "simulation_year",
    "n_sample",
    "model_x",
    "model_y",
    "radius_model_units",
    "inside_model"
)

missing_cols <- setdiff(required_target_cols, names(targets))

if (length(missing_cols) > 0) {
    stop(
        "validation_targets.tsv is missing columns: ",
        paste(missing_cols, collapse = ", ")
    )
}

# Only use targets that lie within the model extent.
targets <- targets[targets$inside_model %in% c(TRUE, "TRUE", "T", 1), ]

cat("Validation targets loaded:", nrow(targets), "\n\n")

# ------------------------------------------------------------
# Find position files
# ------------------------------------------------------------

files <- list.files(
    "sampling",
    pattern = "_positions\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
)

if (length(files) == 0)
    stop("No *_positions.tsv files found under sampling/")

cat("Position files found:", length(files), "\n\n")

# ------------------------------------------------------------
# Parse path / filename
# Example:
# Kexp1p3_D120_rep001_year1951_positions.tsv
# ------------------------------------------------------------

parse_position_file <- function(f) {

    b <- basename(f)

    m <- regexec(
        "^(.*)_(rep[0-9]{3})_year([0-9]{4})_positions\\.tsv$",
        b
    )

    parts <- regmatches(b, m)[[1]]

    if (length(parts) != 4)
        stop("Could not parse filename: ", f)

    if (grepl("/skane_only/", f, fixed = TRUE)) {
        landscape <- "skane_only"
    } else if (grepl("/no_oland/", f, fixed = TRUE)) {
        landscape <- "no_oland"
    } else if (grepl("/full/", f, fixed = TRUE)) {
        landscape <- "full"
    } else {
        landscape <- "unknown"
    }

    list(
        landscape = landscape,
        model      = parts[2],
        rep        = parts[3],
        year       = as.integer(parts[4])
    )
}

# ------------------------------------------------------------
# Calculate sampling pools
# ------------------------------------------------------------

results <- list()
counter <- 1

for (f in files) {

    info <- parse_position_file(f)

    if (info$landscape == "unknown") {
        warning("Skipping file with unknown landscape: ", f)
        next
    }

    # Targets corresponding to this simulation snapshot year.
    tt <- targets[
        targets$simulation_year == info$year,
        ,
        drop = FALSE
    ]

    if (nrow(tt) == 0)
        next

    # Småland is outside the Skåne-only model.
    if (info$landscape == "skane_only") {
        tt <- tt[
            tt$region != "NSmaland",
            ,
            drop = FALSE
        ]
    }

    if (nrow(tt) == 0)
        next

    pos <- read.delim(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    required_pos_cols <- c(
        "model",
        "rep",
        "year",
        "ind_index",
        "sex",
        "age",
        "x",
        "y"
    )

    missing_pos <- setdiff(required_pos_cols, names(pos))

    if (length(missing_pos) > 0) {
        stop(
            "Position file missing columns: ",
            f,
            "\n",
            paste(missing_pos, collapse = ", ")
        )
    }

    for (i in seq_len(nrow(tt))) {

        target <- tt[i, ]

        d <- sqrt(
            (pos$x - target$model_x)^2 +
            (pos$y - target$model_y)^2
        )

        radius <- target$radius_model_units
        inside <- d <= radius

        available <- sum(inside)

        nearest_dist <- if (length(d) > 0) {
            min(d)
        } else {
            NA_real_
        }

        n_female <- sum(pos$sex[inside] == "F")
        n_male   <- sum(pos$sex[inside] == "M")

        results[[counter]] <- data.frame(
            landscape         = info$landscape,
            model             = info$model,
            rep               = info$rep,
            simulation_year   = info$year,

            target             = target$target,
            target_region      = target$region,
            locality           = target$locality,
            empirical_population = target$empirical_population,
            empirical_year     = target$empirical_year,
            match_type         = target$match_type,

            model_x            = target$model_x,
            model_y            = target$model_y,
            radius             = radius,

            total_snapshot_N   = nrow(pos),
            available          = available,
            required           = target$n_sample,
            enough             = available >= target$n_sample,

            nearest_distance   = nearest_dist,

            females_in_pool    = n_female,
            males_in_pool      = n_male,

            source_file        = f,

            stringsAsFactors = FALSE
        )

        counter <- counter + 1
    }
}

if (length(results) == 0)
    stop("No matching target/simulation combinations were found.")

pool <- do.call(rbind, results)

pool <- pool[
    order(
        pool$landscape,
        pool$model,
        pool$target,
        pool$rep
    ),
]

# ------------------------------------------------------------
# Save detailed result
# ------------------------------------------------------------

write.table(
    pool,
    "results/sampling_pool_counts.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------
# Summarise over replicates
# ------------------------------------------------------------

group_key <- interaction(
    pool$landscape,
    pool$model,
    pool$target,
    drop = TRUE
)

groups <- split(pool, group_key)

summary_list <- lapply(groups, function(d) {

    reps_checked <- length(unique(d$rep))

    data.frame(
        landscape          = d$landscape[1],
        model              = d$model[1],
        target             = d$target[1],
        target_region      = d$target_region[1],
        locality           = d$locality[1],
        empirical_year     = d$empirical_year[1],
        simulation_year    = d$simulation_year[1],
        required           = d$required[1],

        reps_checked       = reps_checked,
        complete_5_reps    = reps_checked == 5,

        min_available      = min(d$available),
        median_available   = median(d$available),
        max_available      = max(d$available),

        reps_enough        = sum(d$enough),
        all_checked_enough = all(d$enough),

        all_5_reps_enough  = (
            reps_checked == 5 &&
            all(d$enough)
        ),

        max_nearest_distance = max(d$nearest_distance),

        stringsAsFactors = FALSE
    )
})

summary_out <- do.call(rbind, summary_list)

summary_out <- summary_out[
    order(
        summary_out$landscape,
        summary_out$model,
        summary_out$target
    ),
]

write.table(
    summary_out,
    "results/sampling_pool_summary.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLING POOL SUMMARY\n")
cat("============================================================\n\n")

print(
    summary_out[
        ,
        c(
            "landscape",
            "model",
            "target",
            "required",
            "reps_checked",
            "min_available",
            "median_available",
            "max_available",
            "reps_enough",
            "all_checked_enough"
        )
    ],
    row.names = FALSE
)

cat("\n============================================================\n")
cat("INSUFFICIENT POOLS AMONG CHECKED REPLICATES\n")
cat("============================================================\n\n")

bad <- summary_out[
    !summary_out$all_checked_enough,
    ,
    drop = FALSE
]

if (nrow(bad) == 0) {

    cat("NONE\n")

} else {

    print(
        bad[
            ,
            c(
                "landscape",
                "model",
                "target",
                "required",
                "reps_checked",
                "min_available",
                "reps_enough"
            )
        ],
        row.names = FALSE
    )
}

cat("\nDetailed output:\n")
cat("  results/sampling_pool_counts.tsv\n")
cat("Summary output:\n")
cat("  results/sampling_pool_summary.tsv\n")
