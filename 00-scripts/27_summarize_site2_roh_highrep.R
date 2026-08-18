options(stringsAsFactors = FALSE)

BASE <- path.expand(
    "~/tabea_work/08_population_analysis_final/skane_only/results/regional_roh_site2"
)

OUT <- file.path(
    BASE,
    "combined"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# HELPERS
# ============================================================

mean_na <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) == 0)
        return(NA_real_)

    mean(x)
}

sd_na <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) < 2)
        return(NA_real_)

    sd(x)
}

se_na <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) < 2)
        return(NA_real_)

    sd(x) / sqrt(length(x))
}

# ============================================================
# 1. FIND INPUT FILES
# ============================================================

individual_files <- list.files(
    file.path(BASE, "individuals"),
    pattern = "_individual_ROH_metrics\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
)

regional_files <- list.files(
    file.path(BASE, "summaries"),
    pattern = "_regional_ROH_summary\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
)

segment_files <- list.files(
    file.path(BASE, "segments"),
    pattern = "_ROH_segments_gt100kb\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
)

cat("Individual files :", length(individual_files), "\n")
cat("Regional files   :", length(regional_files), "\n")
cat("Segment files    :", length(segment_files), "\n")

if (
    length(individual_files) == 0 ||
    length(regional_files) == 0 ||
    length(segment_files) == 0
) {
    stop("No ROH result files found.")
}

if (
    length(individual_files) != length(regional_files) ||
    length(individual_files) != length(segment_files)
) {
    stop(
        "Different numbers of individual, regional and segment files."
    )
}

# ============================================================
# 2. READ EVERYTHING
# ============================================================

individual <- do.call(
    rbind,
    lapply(
        individual_files,
        read.delim,
        check.names = FALSE
    )
)

regional <- do.call(
    rbind,
    lapply(
        regional_files,
        read.delim,
        check.names = FALSE
    )
)

segments <- do.call(
    rbind,
    lapply(
        segment_files,
        read.delim,
        check.names = FALSE
    )
)

rownames(individual) <- NULL
rownames(regional) <- NULL
rownames(segments) <- NULL

cat("\nCombined individual rows :", nrow(individual), "\n")
cat("Combined regional rows   :", nrow(regional), "\n")
cat("Combined ROH segments    :", nrow(segments), "\n")

# ============================================================
# 3. DEFINE STATE LABEL
# ============================================================

make_state <- function(scenario, year) {

    if (scenario %in% c("shared_history", "historical")) {

        return(
            paste0(year)
        )
    }

    if (scenario == "status_quo")
        return("2140 SQ")

    if (scenario == "restore_2km")
        return("2140 R2")

    if (scenario == "restore_4km")
        return("2140 R4")

    if (scenario == "restore_6km")
        return("2140 R6")

    paste(
        scenario,
        year
    )
}

individual$state <- mapply(
    make_state,
    individual$scenario,
    individual$year
)

regional$state <- mapply(
    make_state,
    regional$scenario,
    regional$year
)

if (nrow(segments) > 0) {

    segments$state <- mapply(
        make_state,
        segments$scenario,
        segments$year
    )
}

state_order <- c(
    "1900",
    "2020",
    "2140 SQ",
    "2140 R2",
    "2140 R4",
    "2140 R6"
)

region_order <- c(
    "NW",
    "NE",
    "W",
    "E",
    "SE"
)

individual$state <- factor(
    individual$state,
    levels = state_order
)

regional$state <- factor(
    regional$state,
    levels = state_order
)

individual$region <- factor(
    individual$region,
    levels = region_order
)

regional$region <- factor(
    regional$region,
    levels = region_order
)

if (nrow(segments) > 0) {

    segments$state <- factor(
        segments$state,
        levels = state_order
    )

    segments$region <- factor(
        segments$region,
        levels = region_order
    )
}

# ============================================================
# 4. WRITE FULL COMBINED TABLES
# ============================================================

write.table(
    individual,
    file.path(
        OUT,
        "regional_ROH_all_individuals.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    segments,
    file.path(
        OUT,
        "regional_ROH_all_segments_gt100kb.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    regional,
    file.path(
        OUT,
        "regional_ROH_by_replicate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. REGIONAL MEANS ACROSS THE FIVE TRUE REPLICATES
#
# Each row of 'regional' is already:
# one region x one simulation replicate x one state.
#
# Therefore replicate is the unit used here for uncertainty.
# ============================================================

metric_cols <- c(
    "mean_FROH_100kb",
    "mean_n_roh",
    "mean_roh_length_bp",
    "mean_individual_max_roh_bp",
    "max_roh_bp_observed",
    "mean_FROH_100_250kb",
    "mean_FROH_250_500kb",
    "mean_FROH_500kb_1Mb",
    "mean_FROH_gt1Mb"
)

groups <- split(
    regional,
    interaction(
        regional$state,
        regional$region,
        drop = TRUE
    )
)

regional_summary_list <- lapply(
    groups,
    function(g) {

        out <- data.frame(
            state = as.character(g$state[1]),
            year = g$year[1],
            scenario = g$scenario[1],
            region = as.character(g$region[1]),
            n_reps = nrow(g),
            stringsAsFactors = FALSE
        )

        for (v in metric_cols) {

            out[[paste0(v, "_mean")]] <-
                mean_na(g[[v]])

            out[[paste0(v, "_sd")]] <-
                sd_na(g[[v]])

            out[[paste0(v, "_se")]] <-
                se_na(g[[v]])
        }

        out
    }
)

regional_summary <- do.call(
    rbind,
    regional_summary_list
)

regional_summary$state <- factor(
    regional_summary$state,
    levels = state_order
)

regional_summary$region <- factor(
    regional_summary$region,
    levels = region_order
)

regional_summary <- regional_summary[
    order(
        regional_summary$state,
        regional_summary$region
    ),
]

write.table(
    regional_summary,
    file.path(
        OUT,
        "regional_ROH_summary_across_replicates.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. WHOLE-SKANE VALUE PER REPLICATE
#
# Equal-weight mean across the five regional means.
#
# This is an equal-weight mean across the five regional
# site-2 samples. Regional sample sizes may differ slightly
# where fewer than 20 individuals were available.
# ============================================================

rep_groups <- split(
    regional,
    interaction(
        regional$state,
        regional$rep,
        drop = TRUE
    )
)

skane_rep_list <- lapply(
    rep_groups,
    function(g) {

        out <- data.frame(
            state = as.character(g$state[1]),
            year = g$year[1],
            scenario = g$scenario[1],
            rep = g$rep[1],
            n_regions = nrow(g),
            stringsAsFactors = FALSE
        )

        if (nrow(g) != 5)
            stop(
                "Expected five regions for ",
                as.character(g$state[1]),
                " ",
                g$rep[1]
            )

        for (v in metric_cols) {

            out[[v]] <-
                mean_na(
                    g[[v]]
                )
        }

        out
    }
)

skane_rep <- do.call(
    rbind,
    skane_rep_list
)

skane_rep$state <- factor(
    skane_rep$state,
    levels = state_order
)

skane_rep <- skane_rep[
    order(
        skane_rep$state,
        skane_rep$rep
    ),
]

write.table(
    skane_rep,
    file.path(
        OUT,
        "skane_ROH_by_replicate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. WHOLE-SKANE SUMMARY ACROSS REPLICATES
# ============================================================

state_groups <- split(
    skane_rep,
    skane_rep$state,
    drop = TRUE
)

skane_summary_list <- lapply(
    state_groups,
    function(g) {

        out <- data.frame(
            state = as.character(g$state[1]),
            year = g$year[1],
            scenario = g$scenario[1],
            n_reps = nrow(g),
            stringsAsFactors = FALSE
        )

        for (v in metric_cols) {

            out[[paste0(v, "_mean")]] <-
                mean_na(g[[v]])

            out[[paste0(v, "_sd")]] <-
                sd_na(g[[v]])

            out[[paste0(v, "_se")]] <-
                se_na(g[[v]])
        }

        out
    }
)

skane_summary <- do.call(
    rbind,
    skane_summary_list
)

skane_summary$state <- factor(
    skane_summary$state,
    levels = state_order
)

skane_summary <- skane_summary[
    order(skane_summary$state),
]

write.table(
    skane_summary,
    file.path(
        OUT,
        "skane_ROH_summary_across_replicates.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. LENGTH-CLASS CONTRIBUTION
#
# These values quantify how much of mean FROH is contributed
# by each ROH size class.
# ============================================================

length_class <- rbind(

    data.frame(
        regional[, c(
            "state",
            "scenario",
            "rep",
            "year",
            "region"
        )],
        roh_class = "100-250kb",
        FROH_contribution =
            regional$mean_FROH_100_250kb
    ),

    data.frame(
        regional[, c(
            "state",
            "scenario",
            "rep",
            "year",
            "region"
        )],
        roh_class = "250-500kb",
        FROH_contribution =
            regional$mean_FROH_250_500kb
    ),

    data.frame(
        regional[, c(
            "state",
            "scenario",
            "rep",
            "year",
            "region"
        )],
        roh_class = "500kb-1Mb",
        FROH_contribution =
            regional$mean_FROH_500kb_1Mb
    ),

    data.frame(
        regional[, c(
            "state",
            "scenario",
            "rep",
            "year",
            "region"
        )],
        roh_class = ">1Mb",
        FROH_contribution =
            regional$mean_FROH_gt1Mb
    )
)

length_class$roh_class <- factor(
    length_class$roh_class,
    levels = c(
        "100-250kb",
        "250-500kb",
        "500kb-1Mb",
        ">1Mb"
    )
)

write.table(
    length_class,
    file.path(
        OUT,
        "regional_ROH_length_class_by_replicate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 9. DESCRIPTIVE SEGMENT-LENGTH SUMMARY
#
# Note:
# This summarizes actual ROH segments and is therefore
# descriptive. Individuals with more ROHs contribute more
# segments. Inferential summaries remain individual/replicate
# based above.
# ============================================================

if (nrow(segments) > 0) {

    seg_groups <- split(
        segments,
        interaction(
            segments$state,
            segments$rep,
            segments$region,
            drop = TRUE
        )
    )

    seg_summary_list <- lapply(
        seg_groups,
        function(g) {

            data.frame(
                state = as.character(g$state[1]),
                scenario = g$scenario[1],
                rep = g$rep[1],
                year = g$year[1],
                region = as.character(g$region[1]),

                n_segments = nrow(g),

                median_length_bp =
                    median(g$length_bp),

                mean_length_bp =
                    mean(g$length_bp),

                p90_length_bp =
                    unname(
                        quantile(
                            g$length_bp,
                            0.90
                        )
                    ),

                max_length_bp =
                    max(g$length_bp),

                stringsAsFactors = FALSE
            )
        }
    )

    seg_summary <- do.call(
        rbind,
        seg_summary_list
    )

    write.table(
        seg_summary,
        file.path(
            OUT,
            "regional_ROH_segment_length_descriptive.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# ============================================================
# 10. PRINT USEFUL FIRST SUMMARY
# ============================================================

cat("\n============================================\n")
cat("WHOLE-SKANE FROH SUMMARY\n")
cat("============================================\n\n")

print(
    skane_summary[, c(
        "state",
        "n_reps",
        "mean_FROH_100kb_mean",
        "mean_FROH_100kb_sd",
        "mean_n_roh_mean",
        "mean_roh_length_bp_mean",
        "mean_FROH_gt1Mb_mean"
    )],
    row.names = FALSE
)

cat("\n============================================\n")
cat("REGIONAL FROH SUMMARY\n")
cat("============================================\n\n")

print(
    regional_summary[, c(
        "state",
        "region",
        "mean_FROH_100kb_mean",
        "mean_FROH_100kb_sd",
        "mean_n_roh_mean",
        "mean_roh_length_bp_mean",
        "mean_FROH_gt1Mb_mean"
    )],
    row.names = FALSE
)

cat("\nOutputs written to:\n")
cat(OUT, "\n")
cat("\nDONE\n")
