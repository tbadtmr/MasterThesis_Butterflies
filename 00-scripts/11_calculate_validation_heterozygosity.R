options(stringsAsFactors = FALSE)

CHROM_LEN <- 19361589L

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1)
    stop("Usage: Rscript 11_calculate_validation_heterozygosity.R TASK_ID")

task_id <- as.integer(args[1])

# ============================================================
# 1. TASK INFORMATION
# ============================================================

tasks <- read.delim(
    "results/genotype_export_tasks.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

task <- tasks[
    tasks$task_id == task_id,
    ,
    drop = FALSE
]

if (nrow(task) != 1)
    stop("Could not uniquely identify task_id = ", task_id)

landscape <- task$landscape
model     <- task$model
rep       <- task$rep
year      <- task$year

vcf_file <- file.path(
    "validation_genotypes",
    landscape,
    model,
    rep,
    sprintf(
        "%s_%s_%s_year%d.vcf",
        landscape,
        model,
        rep,
        year
    )
)

if (!file.exists(vcf_file))
    stop("VCF not found: ", vcf_file)

# ============================================================
# 2. RELEVANT SAMPLING MANIFEST
# ============================================================

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

cat("============================================\n")
cat("VALIDATION HETEROZYGOSITY\n")
cat("============================================\n")
cat("Task       :", task_id, "\n")
cat("Landscape  :", landscape, "\n")
cat("Model      :", model, "\n")
cat("Rep        :", rep, "\n")
cat("Year       :", year, "\n")
cat("VCF        :", vcf_file, "\n")
cat("Chr length :", CHROM_LEN, "\n")
cat("Manifest rows:", nrow(m), "\n\n")

# ============================================================
# 3. LOCATE VCF HEADER
# ============================================================

con <- file(vcf_file, open = "r")
skip <- 0L

repeat {

    line <- readLines(
        con,
        n = 1,
        warn = FALSE
    )

    if (length(line) == 0) {
        close(con)
        stop("VCF #CHROM header not found.")
    }

    if (startsWith(line, "#CHROM"))
        break

    skip <- skip + 1L
}

close(con)

# ============================================================
# 4. READ VCF
# ============================================================

vcf <- read.delim(
    vcf_file,
    skip = skip,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (ncol(vcf) < 10)
    stop("VCF contains no sample columns.")

if (any(duplicated(vcf$POS)))
    stop("Duplicate genomic positions detected.")

if (any(grepl(",", vcf$ALT, fixed = TRUE)))
    stop("Multiallelic records detected.")

cat("VCF variant records:", nrow(vcf), "\n")
cat("VCF individuals    :", ncol(vcf) - 9L, "\n\n")

# ============================================================
# 5. EXTRACT GENOTYPES
# ============================================================

gt_raw <- as.matrix(
    vcf[, 10:ncol(vcf), drop = FALSE]
)

# Keep only GT field if FORMAT contains other fields.
gt <- sub(
    ":.*$",
    "",
    gt_raw
)

dim(gt) <- dim(gt_raw)

colnames(gt) <- colnames(gt_raw)

allowed <- c(
    "0|0",
    "0|1",
    "1|0",
    "1|1",
    "0/0",
    "0/1",
    "1/0",
    "1/1"
)

bad <- unique(
    gt[!(gt %in% allowed)]
)

if (length(bad) > 0) {
    stop(
        "Unexpected or missing genotype(s): ",
        paste(head(bad, 20), collapse = ", ")
    )
}

# ============================================================
# 6. INDIVIDUAL HETEROZYGOSITY
#
# Nolen et al.:
#
#   heterozygous called genotypes
#   -----------------------------
#      total called genotypes
#
# reported per 1000 bp/genotypes.
#
# Simulated genotypes are complete, so every genomic position
# is callable. Positions absent from the VCF are invariant and
# therefore homozygous.
# ============================================================

het <- (
    gt == "0|1" |
    gt == "1|0" |
    gt == "0/1" |
    gt == "1/0"
)

het_sites <- colSums(het)

H_per_bp <- het_sites / CHROM_LEN

H_per_1000bp <- H_per_bp * 1000

sample_names <- colnames(gt)

if (!all(grepl("^p1:i[0-9]+$", sample_names))) {
    stop(
        "Unexpected VCF sample name format. ",
        "Expected names like p1:i1025."
    )
}

ind_index <- as.integer(
    sub("^p1:i", "", sample_names)
)

individual_H <- data.frame(
    task_id = task_id,
    landscape = landscape,
    model = model,
    rep = rep,
    simulation_year = year,
    sample_name = sample_names,
    ind_index = ind_index,
    heterozygous_sites = as.integer(het_sites),
    chromosome_length = CHROM_LEN,
    H_per_bp = H_per_bp,
    H_per_1000bp = H_per_1000bp,
    stringsAsFactors = FALSE
)

# ============================================================
# 7. VERIFY MANIFEST -> VCF MAPPING
# ============================================================

required_samples <- unique(
    paste0("p1:i", m$ind_index)
)

missing_samples <- setdiff(
    required_samples,
    individual_H$sample_name
)

if (length(missing_samples) > 0) {
    stop(
        "Manifest individuals missing from VCF: ",
        paste(head(missing_samples, 20), collapse = ", ")
    )
}

# ============================================================
# 8. MEAN H FOR EACH EMPIRICAL-MATCHED DRAW
#
# H itself is an individual statistic.
#
# For comparison with empirical population means, calculate
# mean H across the n sampled individuals in each draw.
# ============================================================

group_key <- paste(
    m$target,
    m$draw,
    sep = "___"
)

groups <- split(
    m,
    group_key
)

draw_res <- vector(
    "list",
    length(groups)
)

j <- 1L

for (g in groups) {

    if (anyDuplicated(g$ind_index))
        stop("Duplicate individual detected within draw.")

    sample_names_draw <- paste0(
        "p1:i",
        g$ind_index
    )

    ii <- match(
        sample_names_draw,
        individual_H$sample_name
    )

    if (anyNA(ii))
        stop("Could not map draw individuals to H table.")

    H_values <- individual_H$H_per_1000bp[ii]

    n_ind <- length(H_values)

    if (
        length(unique(g$n_required)) != 1 ||
        n_ind != unique(g$n_required)
    ) {
        stop(
            "Draw size does not match required empirical sample size."
        )
    }

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

        n_individuals = n_ind,

        mean_H_per_1000bp =
            mean(H_values),

        median_H_per_1000bp =
            median(H_values),

        sd_H_among_individuals =
            if (n_ind > 1) sd(H_values) else NA,

        min_H_per_1000bp =
            min(H_values),

        max_H_per_1000bp =
            max(H_values),

        stringsAsFactors = FALSE
    )

    j <- j + 1L
}

draw_H <- do.call(
    rbind,
    draw_res
)

draw_H <- draw_H[
    order(
        draw_H$target,
        draw_H$draw
    ),
    ,
    drop = FALSE
]

# ============================================================
# 9. WRITE OUTPUTS
# ============================================================

dir.create(
    "results/heterozygosity/per_snapshot",
    recursive = TRUE,
    showWarnings = FALSE
)

individual_file <- file.path(
    "results/heterozygosity/per_snapshot",
    sprintf(
        "H_task%03d_individuals.tsv",
        task_id
    )
)

draw_file <- file.path(
    "results/heterozygosity/per_snapshot",
    sprintf(
        "H_task%03d_draws.tsv",
        task_id
    )
)

write.table(
    individual_H,
    individual_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    draw_H,
    draw_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 10. CONSOLE QC
# ============================================================

cat("Individual H estimates :", nrow(individual_H), "\n")
cat("Draw-level estimates   :", nrow(draw_H), "\n")

cat(
    "Mean individual H / 1000bp :",
    mean(individual_H$H_per_1000bp),
    "\n"
)

cat(
    "Mean draw H / 1000bp       :",
    mean(draw_H$mean_H_per_1000bp),
    "\n"
)

cat(
    "Min draw H / 1000bp        :",
    min(draw_H$mean_H_per_1000bp),
    "\n"
)

cat(
    "Max draw H / 1000bp        :",
    max(draw_H$mean_H_per_1000bp),
    "\n"
)

cat("\nWritten:\n")
cat(individual_file, "\n")
cat(draw_file, "\n")
