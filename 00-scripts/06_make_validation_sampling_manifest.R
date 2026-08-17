options(stringsAsFactors = FALSE)

N_DRAWS <- 100
MASTER_SEED <- 20260810

dir.create("results", showWarnings = FALSE, recursive = TRUE)

pool <- read.delim(
    "results/sampling_pool_counts.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Stable ordering = stable assignment of per-group random seeds.
pool <- pool[
    order(
        pool$landscape,
        pool$model,
        pool$target,
        pool$rep
    ),
]

sampled <- list()
failed  <- list()

s <- 1
f <- 1

for (i in seq_len(nrow(pool))) {

    g <- pool[i, ]

    # Keep failed empirical matches explicitly rather than moving the target.
    if (!isTRUE(g$enough)) {

        failed[[f]] <- data.frame(
            landscape       = g$landscape,
            model           = g$model,
            rep             = g$rep,
            target          = g$target,
            simulation_year = g$simulation_year,
            available       = g$available,
            required        = g$required,
            reason          = "insufficient_individuals_within_radius",
            stringsAsFactors = FALSE
        )

        f <- f + 1
        next
    }

    pos <- read.delim(
        g$source_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    d <- sqrt(
        (pos$x - g$model_x)^2 +
        (pos$y - g$model_y)^2
    )

    candidates <- pos[d <= g$radius, , drop = FALSE]
    candidates$distance_to_target <- d[d <= g$radius]

    if (nrow(candidates) != g$available) {
        stop(
            "Candidate count mismatch for ",
            g$landscape, " / ",
            g$model, " / ",
            g$rep, " / ",
            g$target
        )
    }

    # Unique deterministic seed for this model/rep/target combination.
    group_seed <- MASTER_SEED + i
    set.seed(group_seed)

    for (draw in seq_len(N_DRAWS)) {

        chosen <- sample(
            seq_len(nrow(candidates)),
            size = g$required,
            replace = FALSE
        )

        x <- candidates[chosen, , drop = FALSE]

        sampled[[s]] <- data.frame(
            landscape       = g$landscape,
            model           = g$model,
            rep             = g$rep,
            simulation_year = g$simulation_year,
            target          = g$target,
            empirical_year  = g$empirical_year,
            locality        = g$locality,
            draw             = sprintf("draw%03d", draw),
            sampling_seed   = group_seed,
            n_required      = g$required,
            pool_size       = nrow(candidates),
            ind_index       = x$ind_index,
            sex             = x$sex,
            age             = x$age,
            x               = x$x,
            y               = x$y,
            distance        = x$distance_to_target,
            source_file     = g$source_file,
            stringsAsFactors = FALSE
        )

        s <- s + 1
    }
}

manifest <- do.call(rbind, sampled)

write.table(
    manifest,
    "results/validation_sampling_manifest.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

if (length(failed) > 0) {

    failed_out <- do.call(rbind, failed)

    write.table(
        failed_out,
        "results/validation_sampling_unavailable.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# One line per sample draw for convenient QC.
draw_summary <- aggregate(
    ind_index ~ landscape + model + rep +
        simulation_year + target + draw,
    data = manifest,
    FUN = length
)

names(draw_summary)[names(draw_summary) == "ind_index"] <- "n_sampled"

write.table(
    draw_summary,
    "results/validation_sampling_draw_summary.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("============================================\n")
cat("VALIDATION SAMPLING MANIFEST\n")
cat("============================================\n")
cat("Random draws per available target/rep:", N_DRAWS, "\n")
cat("Master seed:", MASTER_SEED, "\n")
cat("Manifest rows:", nrow(manifest), "\n")
cat("Sample draws:", nrow(draw_summary), "\n")
cat("Unavailable target/rep combinations:", length(failed), "\n\n")

cat("Files written:\n")
cat(" results/validation_sampling_manifest.tsv\n")
cat(" results/validation_sampling_draw_summary.tsv\n")

if (length(failed) > 0)
    cat(" results/validation_sampling_unavailable.tsv\n")
