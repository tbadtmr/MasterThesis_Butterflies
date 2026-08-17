options(stringsAsFactors = FALSE)

INDIR  <- "results/heterozygosity/per_snapshot"
OUTDIR <- "results/heterozygosity"

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 1. READ ALL DRAW-LEVEL FILES
# ============================================================

files <- list.files(
    INDIR,
    pattern = "^H_task[0-9]+_draws\\.tsv$",
    full.names = TRUE
)

cat("Draw-level H files found:", length(files), "\n")

if (length(files) != 94)
    stop(
        "Expected 94 draw-level H files, found ",
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
        "Expected 15500 draw-level rows, found ",
        nrow(dat)
    )

# ============================================================
# 2. EMPIRICAL H VALUES
#
# Current Nolen et al. Supplement Table S2.
#
# Units:
# heterozygous sites / 1000 bp
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

    empirical_H = c(
        2.2752,
        2.2240,
        1.9270,
        2.0940,
        2.022333333333333,
        2.2490,
        2.1470
    ),

    stringsAsFactors = FALSE
)

dat$empirical_H <- empirical$empirical_H[
    match(
        dat$target,
        empirical$target
    )
]

dat$sim_empirical_ratio <- (
    dat$mean_H_per_1000bp /
    dat$empirical_H
)

write.table(
    dat,
    file.path(
        OUTDIR,
        "validation_H_all_draws.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 3. SUMMARIZE 100 DRAWS WITHIN EACH REPLICATE x TARGET
#
# The 100 draws represent sampling uncertainty.
# They are NOT independent simulation replicates.
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

groups <- split(
    dat,
    key
)

rep_target <- do.call(
    rbind,
    lapply(
        groups,
        function(g) {

            x <- g$mean_H_per_1000bp

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

                mean_H =
                    mean(x),

                median_H =
                    median(x),

                sd_sampling =
                    sd(x),

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

                min_H =
                    min(x),

                max_H =
                    max(x),

                empirical_H =
                    g$empirical_H[1],

                sim_empirical_ratio =
                    mean(x) /
                    g$empirical_H[1],

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
        "H_by_replicate_target.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 4. MODEL-LEVEL TARGET SUMMARY
#
# First summarize the 100 sampling draws within each actual
# SLiM replicate. Then compare the independent replicates.
# ============================================================

key2 <- interaction(
    rep_target$landscape,
    rep_target$model,
    rep_target$target,
    drop = TRUE,
    lex.order = TRUE
)

groups2 <- split(
    rep_target,
    key2
)

model_target <- do.call(
    rbind,
    lapply(
        groups2,
        function(g) {

            x <- g$mean_H

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                target =
                    g$target[1],

                n_reps =
                    length(x),

                mean_H =
                    mean(x),

                sd_between_reps =
                    if (length(x) > 1)
                        sd(x)
                    else
                        NA,

                min_rep_H =
                    min(x),

                max_rep_H =
                    max(x),

                empirical_H =
                    g$empirical_H[1],

                sim_empirical_ratio =
                    mean(x) /
                    g$empirical_H[1],

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
        "H_by_model_target.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. HISTORICAL -> MODERN RETENTION
#
# Same model replicate only.
#
# As with pi, Smaland is included only where both targets
# actually exist in that replicate.
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

    empirical_hist_H = c(
        2.2490,
        2.2752,
        2.2752,
        2.0940
    ),

    empirical_modern_H = c(
        2.1470,
        2.2240,
        1.9270,
        2.022333333333333
    ),

    empirical_hist_year = c(
        1951,
        1956,
        1956,
        1936
    ),

    empirical_modern_year = c(
        2021,
        2021,
        2022,
        2021
    ),

    stringsAsFactors = FALSE
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

retention_list <- list()

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

        if (
            nrow(h) != 1 ||
            nrow(m) != 1
        )
            next

        retention <- (
            m$mean_H /
            h$mean_H
        )

        sim_years <- (
            m$simulation_year -
            h$simulation_year
        )

        empirical_retention <- (
            cmp$empirical_modern_H /
            cmp$empirical_hist_H
        )

        empirical_years <- (
            cmp$empirical_modern_year -
            cmp$empirical_hist_year
        )

        retention_100 <- (
            retention^(
                100 / sim_years
            )
        )

        empirical_retention_100 <- (
            empirical_retention^(
                100 / empirical_years
            )
        )

        retention_list[[z]] <- data.frame(

            landscape = land,
            model = mod,
            rep = rp,

            comparison =
                cmp$comparison,

            simulation_hist_year =
                h$simulation_year,

            simulation_modern_year =
                m$simulation_year,

            simulation_years =
                sim_years,

            historical_H =
                h$mean_H,

            modern_H =
                m$mean_H,

            retention =
                retention,

            retention_100yr =
                retention_100,

            empirical_hist_H =
                cmp$empirical_hist_H,

            empirical_modern_H =
                cmp$empirical_modern_H,

            empirical_retention =
                empirical_retention,

            empirical_retention_100yr =
                empirical_retention_100,

            error_retention_100yr_pp =
                100 * (
                    retention_100 -
                    empirical_retention_100
                ),

            stringsAsFactors = FALSE
        )

        z <- z + 1L
    }
}

retention_rep <- do.call(
    rbind,
    retention_list
)

write.table(
    retention_rep,
    file.path(
        OUTDIR,
        "H_retention_by_replicate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. RETENTION SUMMARY PER MODEL
# ============================================================

key3 <- interaction(
    retention_rep$landscape,
    retention_rep$model,
    retention_rep$comparison,
    drop = TRUE,
    lex.order = TRUE
)

groups3 <- split(
    retention_rep,
    key3
)

retention_model <- do.call(
    rbind,
    lapply(
        groups3,
        function(g) {

            x <- g$retention
            x100 <- g$retention_100yr

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                comparison =
                    g$comparison[1],

                n_reps =
                    length(x),

                mean_retention =
                    mean(x),

                sd_retention =
                    if (length(x) > 1)
                        sd(x)
                    else
                        NA,

                min_retention =
                    min(x),

                max_retention =
                    max(x),

                mean_retention_100yr =
                    mean(x100),

                sd_retention_100yr =
                    if (length(x100) > 1)
                        sd(x100)
                    else
                        NA,

                empirical_retention =
                    g$empirical_retention[1],

                empirical_retention_100yr =
                    g$empirical_retention_100yr[1],

                error_100yr_percentage_points =
                    100 * (
                        mean(x100) -
                        g$empirical_retention_100yr[1]
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(retention_model) <- NULL

write.table(
    retention_model,
    file.path(
        OUTDIR,
        "H_retention_by_model.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. PROVISIONAL H-ONLY MODEL RANKING
#
# Use only comparisons available across all model families.
# Smaland is excluded from ranking because occupancy is not
# consistently reproduced.
#
# This is NOT final model selection.
# ============================================================

rank_dat <- retention_model[
    retention_model$comparison %in%
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

            errors_pp <- (
                100 * (
                    g$mean_retention_100yr -
                    g$empirical_retention_100yr
                )
            )

            data.frame(
                landscape =
                    g$landscape[1],

                model =
                    g$model[1],

                n_comparisons =
                    nrow(g),

                mean_absolute_error_pp =
                    mean(
                        abs(errors_pp)
                    ),

                RMSE_pp =
                    sqrt(
                        mean(
                            errors_pp^2
                        )
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

ranking <- ranking[
    order(
        ranking$RMSE_pp
    ),
    ,
    drop = FALSE
]

ranking$H_only_rank <- seq_len(
    nrow(ranking)
)

write.table(
    ranking,
    file.path(
        OUTDIR,
        "H_only_model_ranking.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. CONSOLE SUMMARY
# ============================================================

cat("\n============================================\n")
cat("HETEROZYGOSITY SUMMARY COMPLETE\n")
cat("============================================\n")

cat(
    "Draw-level estimates :",
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
    "Retention rep rows   :",
    nrow(retention_rep),
    "\n"
)

cat("\nH-only provisional ranking:\n\n")

print(
    ranking,
    row.names = FALSE
)

cat(
    "\nWritten to:",
    OUTDIR,
    "\n"
)
