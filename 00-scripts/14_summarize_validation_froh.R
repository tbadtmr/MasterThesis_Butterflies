options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
    stop("Usage: Rscript 14_summarize_validation_froh.R TAG")
tag <- args[1]

OUTDIR <- file.path("results", paste0("froh_", tag))
INDIR  <- file.path(OUTDIR, "per_snapshot")

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 1. READ ALL 94 DRAW-LEVEL FROH FILES
# ============================================================

files <- list.files(
    INDIR,
    pattern = "^FROH_task[0-9]+_draws\\.tsv$",
    full.names = TRUE
)

cat("FROH draw files found:", length(files), "\n")

if (length(files) != 94)
    stop(
        "Expected 94 FROH draw files, found ",
        length(files)
    )

dat <- do.call(
    rbind,
    lapply(
        files,
        read.delim,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
)

cat("Total draw-level rows:", nrow(dat), "\n")

if (nrow(dat) != 15500)
    warning(
        "Expected 15500 draw rows, found ",
        nrow(dat)
    )

# ============================================================
# 2. EMPIRICAL FROH
#
# Current Nolen et al. Supplement Table S2.
#
# These are population means calculated from the individual
# FROH >100 kb values.
# ============================================================

empirical <- data.frame(

    target = c(
        "ESkane_hist",
        "ESkane_mod",
        "SESkane_mod",
        "Smaland_hist",
        "Smaland_mod",
        "WSkane_hist",
        "WSkane_mod"
    ),

    empirical_FROH = c(
        0.0084000000,
        0.0161250000,
        0.1537142857,
        0.0060000000,
        0.0685000000,
        0.0082000000,
        0.0613750000
    ),

    stringsAsFactors = FALSE
)

dat$empirical_FROH <- empirical$empirical_FROH[
    match(
        dat$target,
        empirical$target
    )
]

write.table(
    dat,
    file.path(
        OUTDIR,
        "validation_FROH_all_draws.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 3. SUMMARIZE 100 DRAWS WITHIN EACH TRUE SLIM REPLICATE
#
# Important:
# the 100 draws represent sampling uncertainty.
# They are NOT independent model replicates.
# ============================================================

key <- interaction(
    dat$landscape,
    dat$model,
    dat$rep,
    dat$simulation_year,
    dat$target,
    drop = TRUE,
    lex.order = TRUE
)

groups <- split(dat, key)

rep_target <- do.call(
    rbind,
    lapply(
        groups,
        function(g) {

            x <- g$mean_FROH_chr9

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                rep =
                    g$rep[1],

                simulation_year =
                    g$simulation_year[1],

                target =
                    g$target[1],

                empirical_year =
                    g$empirical_year[1],

                locality =
                    g$locality[1],

                n_draws =
                    length(x),

                mean_FROH =
                    mean(x),

                median_draw_FROH =
                    median(x),

                sd_sampling =
                    sd(x),

                q025_sampling =
                    unname(
                        quantile(x, 0.025)
                    ),

                q975_sampling =
                    unname(
                        quantile(x, 0.975)
                    ),

                min_draw_FROH =
                    min(x),

                max_draw_FROH =
                    max(x),

                empirical_FROH =
                    g$empirical_FROH[1],

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(rep_target) <- NULL

write.table(
    rep_target,
    file.path(
        OUTDIR,
        "FROH_by_replicate_target.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 4. MODEL-LEVEL TARGET SUMMARY
#
# First collapse the 100 draws in each actual replicate.
# Then summarize across the independent SLiM replicates.
# ============================================================

key2 <- interaction(
    rep_target$landscape,
    rep_target$model,
    rep_target$target,
    drop = TRUE,
    lex.order = TRUE
)

groups2 <- split(rep_target, key2)

model_target <- do.call(
    rbind,
    lapply(
        groups2,
        function(g) {

            x <- g$mean_FROH

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                target =
                    g$target[1],

                n_reps =
                    length(x),

                mean_FROH =
                    mean(x),

                sd_between_reps =
                    if (length(x) > 1)
                        sd(x)
                    else
                        NA,

                min_rep_FROH =
                    min(x),

                max_rep_FROH =
                    max(x),

                empirical_FROH =
                    g$empirical_FROH[1],

                absolute_error =
                    mean(x) -
                    g$empirical_FROH[1],

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(model_target) <- NULL

write.table(
    model_target,
    file.path(
        OUTDIR,
        "FROH_by_model_target.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. HISTORICAL -> MODERN CHANGE IN FROH
#
# Primary statistic:
#
# delta_FROH = modern - historical
#
# We DO NOT use ratios because historical FROH can be very
# close to zero.
# ============================================================

comparisons <- data.frame(

    comparison = c(
        "W_Skane",
        "E_Skane",
        "Ehist_to_SEmodern",
        "Smaland"
    ),

    hist_target = c(
        "WSkane_hist",
        "ESkane_hist",
        "ESkane_hist",
        "Smaland_hist"
    ),

    modern_target = c(
        "WSkane_mod",
        "ESkane_mod",
        "SESkane_mod",
        "Smaland_mod"
    ),

    empirical_hist_FROH = c(
        0.0082000000,
        0.0084000000,
        0.0084000000,
        0.0060000000
    ),

    empirical_modern_FROH = c(
        0.0613750000,
        0.0161250000,
        0.1537142857,
        0.0685000000
    ),

    stringsAsFactors = FALSE
)

comparisons$empirical_delta_FROH <- (
    comparisons$empirical_modern_FROH -
    comparisons$empirical_hist_FROH
)

ids <- unique(
    rep_target[
        c(
            "landscape",
            "model",
            "rep"
        )
    ]
)

change_list <- list()
z <- 1L

for (i in seq_len(nrow(ids))) {

    land <- ids$landscape[i]
    mod  <- ids$model[i]
    rp   <- ids$rep[i]

    sub <- rep_target[
        rep_target$landscape == land &
        rep_target$model == mod &
        rep_target$rep == rp,
        ,
        drop = FALSE
    ]

    for (j in seq_len(nrow(comparisons))) {

        cmp <- comparisons[j, ]

        h <- sub[
            sub$target ==
                cmp$hist_target,
            ,
            drop = FALSE
        ]

        m <- sub[
            sub$target ==
                cmp$modern_target,
            ,
            drop = FALSE
        ]

        # If either target did not exist in the simulation,
        # no genetic comparison is forced.
        if (
            nrow(h) != 1 ||
            nrow(m) != 1
        )
            next

        delta <- (
            m$mean_FROH -
            h$mean_FROH
        )

        change_list[[z]] <- data.frame(

            landscape = land,
            model = mod,
            rep = rp,

            comparison =
                cmp$comparison,

            simulation_hist_year =
                h$simulation_year,

            simulation_modern_year =
                m$simulation_year,

            historical_FROH =
                h$mean_FROH,

            modern_FROH =
                m$mean_FROH,

            delta_FROH =
                delta,

            empirical_historical_FROH =
                cmp$empirical_hist_FROH,

            empirical_modern_FROH =
                cmp$empirical_modern_FROH,

            empirical_delta_FROH =
                cmp$empirical_delta_FROH,

            delta_error =
                delta -
                cmp$empirical_delta_FROH,

            stringsAsFactors = FALSE
        )

        z <- z + 1L
    }
}

change_rep <- do.call(
    rbind,
    change_list
)

write.table(
    change_rep,
    file.path(
        OUTDIR,
        "FROH_change_by_replicate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. CHANGE SUMMARY PER MODEL
# ============================================================

key3 <- interaction(
    change_rep$landscape,
    change_rep$model,
    change_rep$comparison,
    drop = TRUE,
    lex.order = TRUE
)

groups3 <- split(change_rep, key3)

change_model <- do.call(
    rbind,
    lapply(
        groups3,
        function(g) {

            x <- g$delta_FROH

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                comparison =
                    g$comparison[1],

                n_reps =
                    length(x),

                mean_historical_FROH =
                    mean(g$historical_FROH),

                mean_modern_FROH =
                    mean(g$modern_FROH),

                mean_delta_FROH =
                    mean(x),

                sd_delta_FROH =
                    if (length(x) > 1)
                        sd(x)
                    else
                        NA,

                min_delta_FROH =
                    min(x),

                max_delta_FROH =
                    max(x),

                empirical_historical_FROH =
                    g$empirical_historical_FROH[1],

                empirical_modern_FROH =
                    g$empirical_modern_FROH[1],

                empirical_delta_FROH =
                    g$empirical_delta_FROH[1],

                error_delta_FROH =
                    mean(x) -
                    g$empirical_delta_FROH[1],

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(change_model) <- NULL

write.table(
    change_model,
    file.path(
        OUTDIR,
        "FROH_change_by_model.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. PROVISIONAL FROH-ONLY MODEL RANKING
#
# Use the three Skane comparisons that are available for all
# parameterisations.
#
# Ranking criterion:
# RMSE of absolute error in delta FROH.
#
# This is NOT final model selection.
# ============================================================

rank_dat <- change_model[
    change_model$comparison %in%
        c(
            "W_Skane",
            "E_Skane",
            "Ehist_to_SEmodern"
        ),
    ,
    drop = FALSE
]

rank_key <- interaction(
    rank_dat$landscape,
    rank_dat$model,
    drop = TRUE,
    lex.order = TRUE
)

rank_groups <- split(
    rank_dat,
    rank_key
)

ranking <- do.call(
    rbind,
    lapply(
        rank_groups,
        function(g) {

            errors <- (
                g$mean_delta_FROH -
                g$empirical_delta_FROH
            )

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                n_comparisons =
                    nrow(g),

                mean_absolute_error =
                    mean(abs(errors)),

                RMSE =
                    sqrt(mean(errors^2)),

                stringsAsFactors = FALSE
            )
        }
    )
)

ranking <- ranking[
    order(ranking$RMSE),
    ,
    drop = FALSE
]

ranking$FROH_only_rank <-
    seq_len(nrow(ranking))

write.table(
    ranking,
    file.path(
        OUTDIR,
        "FROH_only_model_ranking.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. CONSOLE SUMMARY
# ============================================================

cat("\n============================================\n")
cat("FROH SUMMARY COMPLETE\n")
cat("============================================\n")

cat(
    "Draw-level rows      :",
    nrow(dat),
    "\n"
)

cat(
    "Replicate-target rows:",
    nrow(rep_target),
    "\n"
)

cat(
    "Model-target rows    :",
    nrow(model_target),
    "\n"
)

cat(
    "Change replicate rows:",
    nrow(change_rep),
    "\n"
)

cat("\nEmpirical delta FROH targets:\n\n")

print(
    comparisons[
        ,
        c(
            "comparison",
            "empirical_hist_FROH",
            "empirical_modern_FROH",
            "empirical_delta_FROH"
        )
    ],
    row.names = FALSE
)

cat("\nFROH-only provisional ranking:\n\n")

print(
    ranking,
    row.names = FALSE
)

cat(
    "\nWritten to:",
    OUTDIR,
    "\n"
)
