options(stringsAsFactors = FALSE)

CHROM_LEN <- 19361589L
MIN_ROH   <- 100000L

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2)
    stop("Usage: Rscript 13_parse_validation_froh.R TASK_ID TAG")

task_id <- as.integer(args[1])

tag <- args[2]
outroot_raw  <- file.path("results", paste0("froh_", tag), "raw")
outroot_snap <- file.path("results", paste0("froh_", tag), "per_snapshot")

tasks <- read.delim(
    "results/genotype_export_tasks.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

task <- tasks[tasks$task_id == task_id, , drop = FALSE]

if (nrow(task) != 1)
    stop("Could not uniquely identify task ", task_id)

landscape <- task$landscape
model     <- task$model
rep       <- task$rep
year      <- task$year

raw_file <- file.path(
    outroot_raw,
    sprintf("roh_task%03d.txt", task_id)
)

if (!file.exists(raw_file))
    stop("ROH file not found: ", raw_file)

H_file <- file.path(
    "results/heterozygosity/per_snapshot",
    sprintf("H_task%03d_individuals.tsv", task_id)
)

if (!file.exists(H_file))
    stop("H file not found: ", H_file)

H <- read.delim(
    H_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# ------------------------------------------------------------
# Read BCFtools region output
# ------------------------------------------------------------

lines <- readLines(raw_file, warn = FALSE)
rg_lines <- lines[grepl("^RG\\t", lines)]

if (length(rg_lines) == 0) {

    roh <- data.frame(
        type = character(),
        sample_name = character(),
        chromosome = character(),
        start = integer(),
        end = integer(),
        length_bp = integer(),
        n_markers = integer(),
        quality = numeric(),
        stringsAsFactors = FALSE
    )

} else {

    roh <- read.delim(
        text = paste(rg_lines, collapse = "\n"),
        header = FALSE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE
    )

    if (ncol(roh) != 8)
        stop("Unexpected number of ROH columns: ", ncol(roh))

    names(roh) <- c(
        "type",
        "sample_name",
        "chromosome",
        "start",
        "end",
        "length_bp",
        "n_markers",
        "quality"
    )
}

cat("============================================\n")
cat("VALIDATION FROH\n")
cat("============================================\n")
cat("Task      :", task_id, "\n")
cat("Landscape :", landscape, "\n")
cat("Model     :", model, "\n")
cat("Rep       :", rep, "\n")
cat("Year      :", year, "\n")
cat("Chr length:", CHROM_LEN, "\n")
cat("ROH cutoff: >", MIN_ROH, "bp\n")
cat("All ROH regions:", nrow(roh), "\n")

# Nolen: ROH >100 kb
roh100 <- roh[
    roh$length_bp > MIN_ROH,
    ,
    drop = FALSE
]

cat("ROH >100 kb:", nrow(roh100), "\n")

# ------------------------------------------------------------
# Individual FROH
# ------------------------------------------------------------

individual <- H[, c(
    "sample_name",
    "ind_index",
    "heterozygous_sites",
    "H_per_1000bp"
)]

individual$n_roh_gt100kb <- 0L
individual$total_roh_bp  <- 0
individual$max_roh_bp    <- 0

if (nrow(roh100) > 0) {

    unknown <- setdiff(
        unique(roh100$sample_name),
        individual$sample_name
    )

    if (length(unknown) > 0)
        stop(
            "Unknown sample names in ROH output: ",
            paste(head(unknown, 20), collapse = ", ")
        )

    by_sample <- split(
        roh100,
        roh100$sample_name
    )

    for (s in names(by_sample)) {

        g <- by_sample[[s]]

        ii <- match(
            s,
            individual$sample_name
        )

        individual$n_roh_gt100kb[ii] <- nrow(g)
        individual$total_roh_bp[ii]  <- sum(g$length_bp)
        individual$max_roh_bp[ii]    <- max(g$length_bp)
    }
}

individual$FROH_chr9 <- (
    individual$total_roh_bp /
    CHROM_LEN
)

if (any(individual$FROH_chr9 > 1 + 1e-8))
    stop("FROH > 1 detected.")

individual$task_id <- task_id
individual$landscape <- landscape
individual$model <- model
individual$rep <- rep
individual$simulation_year <- year

individual <- individual[, c(
    "task_id",
    "landscape",
    "model",
    "rep",
    "simulation_year",
    "sample_name",
    "ind_index",
    "heterozygous_sites",
    "H_per_1000bp",
    "n_roh_gt100kb",
    "total_roh_bp",
    "max_roh_bp",
    "FROH_chr9"
)]

# ------------------------------------------------------------
# Sampling manifest
# ------------------------------------------------------------

manifest <- read.delim(
    "results/validation_sampling_manifest.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

m <- manifest[
    manifest$landscape == landscape &
    manifest$model == model &
    manifest$rep == rep &
    manifest$simulation_year == year,
    ,
    drop = FALSE
]

if (nrow(m) == 0)
    stop("No manifest rows found for task ", task_id)

key <- paste(
    m$target,
    m$draw,
    sep = "___"
)

groups <- split(m, key)

draw_res <- vector(
    "list",
    length(groups)
)

j <- 1L

for (g in groups) {

    if (anyDuplicated(g$ind_index))
        stop("Duplicate individual within draw.")

    ii <- match(
        g$ind_index,
        individual$ind_index
    )

    if (anyNA(ii))
        stop("Could not map draw individuals to FROH table.")

    x <- individual$FROH_chr9[ii]

    if (
        length(unique(g$n_required)) != 1 ||
        length(x) != unique(g$n_required)
    )
        stop("Incorrect draw sample size.")

    draw_res[[j]] <- data.frame(
        task_id = task_id,
        landscape = landscape,
        model = model,
        rep = rep,
        simulation_year = year,
        target = g$target[1],
        empirical_year = g$empirical_year[1],
        locality = g$locality[1],
        draw = g$draw[1],
        n_individuals = length(x),

        mean_FROH_chr9 = mean(x),
        median_FROH_chr9 = median(x),

        sd_FROH_among_individuals =
            if (length(x) > 1)
                sd(x)
            else
                NA,

        min_FROH_chr9 = min(x),
        max_FROH_chr9 = max(x),

        stringsAsFactors = FALSE
    )

    j <- j + 1L
}

draws <- do.call(
    rbind,
    draw_res
)

draws <- draws[
    order(
        draws$target,
        draws$draw
    ),
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

dir.create(
    outroot_snap,
    recursive = TRUE,
    showWarnings = FALSE
)

individual_file <- file.path(
    outroot_snap,
    sprintf("FROH_task%03d_individuals.tsv", task_id)
)

draw_file <- file.path(
    outroot_snap,
    sprintf("FROH_task%03d_draws.tsv", task_id)
)

write.table(
    individual,
    individual_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    draws,
    draw_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

cat("\nIndividuals:", nrow(individual), "\n")
cat("Draws      :", nrow(draws), "\n")

cat("\nIndividual FROH summary:\n")
print(summary(individual$FROH_chr9))

cat(
    "\nFROH = 0  :",
    sum(individual$FROH_chr9 == 0),
    "\n"
)

cat(
    "FROH >0.1:",
    sum(individual$FROH_chr9 > 0.1),
    "\n"
)

cat(
    "FROH >0.5:",
    sum(individual$FROH_chr9 > 0.5),
    "\n"
)

cat(
    "FROH >0.9:",
    sum(individual$FROH_chr9 > 0.9),
    "\n"
)

cat("\nLowest-H individuals and FROH:\n")

tmp <- individual[
    order(individual$H_per_1000bp),
    c(
        "ind_index",
        "H_per_1000bp",
        "n_roh_gt100kb",
        "total_roh_bp",
        "max_roh_bp",
        "FROH_chr9"
    )
]

print(
    head(tmp, 10),
    row.names = FALSE
)

cat("\nWritten:\n")
cat(individual_file, "\n")
cat(draw_file, "\n")
