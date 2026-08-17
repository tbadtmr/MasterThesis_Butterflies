options(stringsAsFactors = FALSE)

manifest <- read.delim(
    "results/validation_sampling_manifest.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

out_root <- "genotype_requests"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Stable ordering
manifest <- manifest[
    order(
        manifest$landscape,
        manifest$model,
        manifest$rep,
        manifest$simulation_year,
        manifest$ind_index
    ),
]

key <- interaction(
    manifest$landscape,
    manifest$model,
    manifest$rep,
    manifest$simulation_year,
    drop = TRUE,
    lex.order = TRUE
)

groups <- split(manifest, key)

summary_list <- list()
k <- 1

for (d in groups) {

    landscape <- d$landscape[1]
    model     <- d$model[1]
    rep       <- d$rep[1]
    year      <- d$simulation_year[1]

    # One row per unique individual needed from this snapshot
    u <- d[
        !duplicated(d$ind_index),
        c(
            "ind_index",
            "sex",
            "age",
            "x",
            "y"
        )
    ]

    u <- u[order(u$ind_index), ]

    out_dir <- file.path(
        out_root,
        landscape,
        model,
        rep
    )

    dir.create(
        out_dir,
        showWarnings = FALSE,
        recursive = TRUE
    )

    out_file <- file.path(
        out_dir,
        sprintf("year%d_indices.tsv", year)
    )

    write.table(
        u,
        out_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    summary_list[[k]] <- data.frame(
        landscape        = landscape,
        model            = model,
        rep              = rep,
        simulation_year  = year,
        targets           = paste(sort(unique(d$target)), collapse = ","),
        manifest_rows     = nrow(d),
        unique_individuals = nrow(u),
        request_file      = out_file,
        stringsAsFactors = FALSE
    )

    k <- k + 1
}

summary_out <- do.call(rbind, summary_list)

summary_out <- summary_out[
    order(
        summary_out$landscape,
        summary_out$model,
        summary_out$rep,
        summary_out$simulation_year
    ),
]

write.table(
    summary_out,
    "results/genotype_request_summary.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("============================================\n")
cat("GENOTYPE REQUESTS\n")
cat("============================================\n")
cat("Snapshot request files:", nrow(summary_out), "\n")
cat("Total unique individual requests across snapshots:",
    sum(summary_out$unique_individuals), "\n")
cat("\nWritten:\n")
cat("  genotype_requests/\n")
cat("  results/genotype_request_summary.tsv\n")
