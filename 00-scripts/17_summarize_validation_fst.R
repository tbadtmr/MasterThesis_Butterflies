options(stringsAsFactors = FALSE)

INDIR  <- "results/fst/per_modelrep"
OUTDIR <- "results/fst"

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 1. READ ALL 30 MODEL-REPLICATE FILES
# ============================================================

files <- list.files(
    INDIR,
    pattern = "^FST_job[0-9]+_draws\\.tsv$",
    full.names = TRUE
)

cat("FST files found:", length(files), "\n")

if (length(files) != 30)
    stop(
        "Expected 30 FST files, found ",
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

rownames(dat) <- NULL

cat("Total draw-level rows:", nrow(dat), "\n")

if (nrow(dat) != 12000)
    stop(
        "Expected 12000 draw rows, found ",
        nrow(dat)
    )

# Check 100 draws in every model-replicate-comparison
check_key <- interaction(
    dat$landscape,
    dat$model,
    dat$rep,
    dat$comparison,
    drop = TRUE
)

draw_counts <- table(check_key)

if (any(draw_counts != 100))
    stop(
        "Not every replicate/comparison contains 100 draws."
    )

write.table(
    dat,
    file.path(
        OUTDIR,
        "validation_FST_all_draws.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 2. COLLAPSE 100 DRAWS WITHIN EACH TRUE SLIM REPLICATE
#
# 100 draws = sampling uncertainty
# 5 SLiM reps = independent simulation replicates
# ============================================================

key <- interaction(
    dat$landscape,
    dat$model,
    dat$rep,
    dat$comparison,
    drop = TRUE,
    lex.order = TRUE
)

groups <- split(dat, key)

rep_comparison <- do.call(
    rbind,
    lapply(
        groups,
        function(g) {

            x <- g$mean_window_FST

            data.frame(
                landscape = g$landscape[1],
                model = g$model[1],
                rep = g$rep[1],
                comparison = g$comparison[1],

                target1 = g$target1[1],
                year1 = g$year1[1],
                n1 = g$n1[1],

                target2 = g$target2[1],
                year2 = g$year2[1],
                n2 = g$n2[1],

                n_draws = length(x),

                mean_FST = mean(x),
                median_FST = median(x),

                sd_sampling = sd(x),

                q025_sampling =
                    unname(
                        quantile(
                            x,
                            0.025
                        )
                    ),

                q975_sampling =
                    unname(
                        quantile(
                            x,
                            0.975
                        )
                    ),

                min_draw_FST = min(x),
                max_draw_FST = max(x),

                mean_component_FST =
                    mean(
                        g$component_FST
                    ),

                mean_n_windows =
                    mean(
                        g$n_windows
                    ),

                mean_n_snps =
                    mean(
                        g$n_snps
                    ),

                empirical_FST =
                    g$empirical_FST[1],

                error =
                    mean(x) -
                    g$empirical_FST[1],

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(rep_comparison) <- NULL

write.table(
    rep_comparison,
    file.path(
        OUTDIR,
        "FST_by_replicate_comparison.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 3. SUMMARIZE ACROSS THE 5 TRUE REPLICATES PER MODEL
# ============================================================

key2 <- interaction(
    rep_comparison$landscape,
    rep_comparison$model,
    rep_comparison$comparison,
    drop = TRUE,
    lex.order = TRUE
)

groups2 <- split(
    rep_comparison,
    key2
)

model_comparison <- do.call(
    rbind,
    lapply(
        groups2,
        function(g) {

            x <- g$mean_FST

            data.frame(
                landscape = g$landscape[1],
                model = g$model[1],
                comparison = g$comparison[1],

                n_reps = length(x),

                mean_FST =
                    mean(x),

                sd_between_reps =
                    if (length(x) > 1)
                        sd(x)
                    else
                        NA,

                min_rep_FST =
                    min(x),

                max_rep_FST =
                    max(x),

                empirical_FST =
                    g$empirical_FST[1],

                absolute_error =
                    mean(x) -
                    g$empirical_FST[1],

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(model_comparison) <- NULL

write.table(
    model_comparison,
    file.path(
        OUTDIR,
        "FST_by_model_comparison.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 4. DIRECT HISTORICAL -> MODERN E-W CHANGE
#
# Empirical:
# historical E-W = 0.02874
# modern E-W     = 0.04940
# delta          = 0.02066
# ============================================================

ids <- unique(
    rep_comparison[
        ,
        c(
            "landscape",
            "model",
            "rep"
        )
    ]
)

delta_list <- list()
z <- 1L

for (i in seq_len(nrow(ids))) {

    sub <- rep_comparison[
        rep_comparison$landscape ==
            ids$landscape[i] &
        rep_comparison$model ==
            ids$model[i] &
        rep_comparison$rep ==
            ids$rep[i],
        ,
        drop = FALSE
    ]

    h <- sub[
        sub$comparison ==
            "historical_EW",
        ,
        drop = FALSE
    ]

    m <- sub[
        sub$comparison ==
            "modern_EW",
        ,
        drop = FALSE
    ]

    if (
        nrow(h) != 1 ||
        nrow(m) != 1
    )
        stop(
            "Missing E-W comparison for ",
            ids$landscape[i], " ",
            ids$model[i], " ",
            ids$rep[i]
        )

    delta <- (
        m$mean_FST -
        h$mean_FST
    )

    delta_list[[z]] <- data.frame(
        landscape =
            ids$landscape[i],

        model =
            ids$model[i],

        rep =
            ids$rep[i],

        historical_EW_FST =
            h$mean_FST,

        modern_EW_FST =
            m$mean_FST,

        delta_FST =
            delta,

        empirical_historical_EW =
            0.02874,

        empirical_modern_EW =
            0.04940,

        empirical_delta_FST =
            0.02066,

        delta_error =
            delta - 0.02066,

        stringsAsFactors = FALSE
    )

    z <- z + 1L
}

delta_rep <- do.call(
    rbind,
    delta_list
)

write.table(
    delta_rep,
    file.path(
        OUTDIR,
        "FST_EW_change_by_replicate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. E-W CHANGE ACROSS REPLICATES
# ============================================================

key3 <- interaction(
    delta_rep$landscape,
    delta_rep$model,
    drop = TRUE,
    lex.order = TRUE
)

groups3 <- split(
    delta_rep,
    key3
)

delta_model <- do.call(
    rbind,
    lapply(
        groups3,
        function(g) {

            x <- g$delta_FST

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                n_reps =
                    nrow(g),

                mean_historical_EW_FST =
                    mean(
                        g$historical_EW_FST
                    ),

                mean_modern_EW_FST =
                    mean(
                        g$modern_EW_FST
                    ),

                mean_delta_FST =
                    mean(x),

                sd_delta_FST =
                    sd(x),

                min_delta_FST =
                    min(x),

                max_delta_FST =
                    max(x),

                empirical_historical_EW =
                    0.02874,

                empirical_modern_EW =
                    0.04940,

                empirical_delta_FST =
                    0.02066,

                error_delta_FST =
                    mean(x) -
                    0.02066,

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(delta_model) <- NULL

write.table(
    delta_model,
    file.path(
        OUTDIR,
        "FST_EW_change_by_model.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. PROVISIONAL FST-ONLY MODEL RANKING
#
# Uses all four directly observed empirical FST benchmarks:
#
# historical E-W
# modern E-W
# modern E-SE
# modern W-SE
#
# This is diagnostic only, not final model selection.
# ============================================================

rank_key <- interaction(
    model_comparison$landscape,
    model_comparison$model,
    drop = TRUE,
    lex.order = TRUE
)

rank_groups <- split(
    model_comparison,
    rank_key
)

ranking <- do.call(
    rbind,
    lapply(
        rank_groups,
        function(g) {

            errors <- (
                g$mean_FST -
                g$empirical_FST
            )

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                n_comparisons =
                    nrow(g),

                mean_absolute_error =
                    mean(
                        abs(errors)
                    ),

                RMSE =
                    sqrt(
                        mean(
                            errors^2
                        )
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

ranking <- ranking[
    order(
        ranking$RMSE
    ),
    ,
    drop = FALSE
]

ranking$FST_only_rank <-
    seq_len(nrow(ranking))

write.table(
    ranking,
    file.path(
        OUTDIR,
        "FST_only_model_ranking.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. QC / CONSOLE SUMMARY
# ============================================================

cat("\n============================================\n")
cat("FST SUMMARY COMPLETE\n")
cat("============================================\n")

cat(
    "Draw-level rows       :",
    nrow(dat),
    "\n"
)

cat(
    "Rep-comparison rows   :",
    nrow(rep_comparison),
    "\n"
)

cat(
    "Model-comparison rows :",
    nrow(model_comparison),
    "\n"
)

cat(
    "E-W delta rep rows    :",
    nrow(delta_rep),
    "\n\n"
)

cat(
    "Model-level FST:\n\n"
)

print(
    model_comparison[
        order(
            model_comparison$landscape,
            model_comparison$model,
            model_comparison$comparison
        ),
    ],
    row.names = FALSE
)

cat(
    "\nE-W temporal change:\n\n"
)

print(
    delta_model,
    row.names = FALSE
)

cat(
    "\nFST-only provisional ranking:\n\n"
)

print(
    ranking,
    row.names = FALSE
)

cat(
    "\nWritten to:",
    OUTDIR,
    "\n"
)
