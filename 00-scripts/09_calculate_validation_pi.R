options(stringsAsFactors = FALSE)

CHROM_LEN  <- 19361589L
WINDOW_SIZE <- 50000L

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1)
    stop("Usage: Rscript 09_calculate_validation_pi.R TASK_ID")

task_id <- as.integer(args[1])

# ============================================================
# TASK
# ============================================================

tasks <- read.delim(
    "results/genotype_export_tasks.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

task <- tasks[tasks$task_id == task_id, , drop = FALSE]

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
# MANIFEST
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
    stop("No manifest rows for task ", task_id)

cat("============================================\n")
cat("PIXY-LIKE VALIDATION PI\n")
cat("============================================\n")
cat("Task       :", task_id, "\n")
cat("Landscape  :", landscape, "\n")
cat("Model      :", model, "\n")
cat("Rep        :", rep, "\n")
cat("Year       :", year, "\n")
cat("VCF        :", vcf_file, "\n")
cat("Window size:", WINDOW_SIZE, "\n")
cat("Chr length :", CHROM_LEN, "\n\n")

# ============================================================
# READ VCF
# ============================================================

con <- file(vcf_file, open = "r")
skip <- 0L

repeat {
    line <- readLines(con, n = 1, warn = FALSE)

    if (length(line) == 0) {
        close(con)
        stop("VCF #CHROM header not found")
    }

    if (startsWith(line, "#CHROM"))
        break

    skip <- skip + 1L
}

close(con)

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

if (any(duplicated(vcf$POS)))
    stop("Duplicate genomic positions detected")

if (any(grepl(",", vcf$ALT, fixed = TRUE)))
    stop("Multiallelic records detected")

cat("VCF variants:", nrow(vcf), "\n")
cat("VCF samples :", ncol(vcf) - 9L, "\n")

# ============================================================
# GENOTYPES -> ALT DOSAGE
# ============================================================

gt_raw <- as.matrix(
    vcf[, 10:ncol(vcf), drop = FALSE]
)

gt <- sub(":.*$", "", gt_raw)
dim(gt) <- dim(gt_raw)

dosage <- matrix(
    NA_integer_,
    nrow = nrow(gt),
    ncol = ncol(gt),
    dimnames = list(NULL, colnames(gt_raw))
)

dosage[gt == "0|0" | gt == "0/0"] <- 0L

dosage[
    gt == "0|1" |
    gt == "1|0" |
    gt == "0/1" |
    gt == "1/0"
] <- 1L

dosage[gt == "1|1" | gt == "1/1"] <- 2L

if (anyNA(dosage)) {
    bad <- unique(gt[is.na(dosage)])

    stop(
        "Unrecognised genotype(s): ",
        paste(head(bad, 20), collapse = ", ")
    )
}

# ============================================================
# 50 KB WINDOWS
#
# VCF positions are 1-based:
#   1-50000
#   50001-100000
#   ...
# ============================================================

n_windows <- ceiling(CHROM_LEN / WINDOW_SIZE)

variant_window <- ((vcf$POS - 1L) %/% WINDOW_SIZE) + 1L

window_start <- ((seq_len(n_windows) - 1L) * WINDOW_SIZE) + 1L
window_end   <- pmin(seq_len(n_windows) * WINDOW_SIZE, CHROM_LEN)
window_length <- window_end - window_start + 1L

cat("Number of 50kb windows:", n_windows, "\n")
cat("Last window length    :", tail(window_length, 1), "\n\n")

# ============================================================
# VERIFY MANIFEST -> VCF
# ============================================================

required_samples <- unique(
    paste0("p1:i", m$ind_index)
)

missing_samples <- setdiff(
    required_samples,
    colnames(dosage)
)

if (length(missing_samples) > 0)
    stop(
        "Manifest individuals missing from VCF: ",
        paste(head(missing_samples, 20), collapse = ", ")
    )

# ============================================================
# CALCULATE PI FOR EACH TARGET x DRAW
# ============================================================

groups <- split(
    m,
    paste(m$target, m$draw, sep = "___")
)

res <- vector("list", length(groups))

j <- 1L

for (g in groups) {

    if (anyDuplicated(g$ind_index))
        stop("Duplicate individual within a draw")

    sample_names <- paste0("p1:i", g$ind_index)

    cols <- match(
        sample_names,
        colnames(dosage)
    )

    if (anyNA(cols))
        stop("Could not map draw individuals to VCF")

    n_ind <- length(cols)

    if (length(unique(g$n_required)) != 1 ||
        n_ind != unique(g$n_required))
        stop("Draw size does not match requested empirical sample size")

    n_chr <- 2L * n_ind

    # Number of ALT chromosomes at each variant position.
    k <- rowSums(
        dosage[, cols, drop = FALSE]
    )

    # Unbiased average pairwise difference at each position.
    #
    # With complete diploid genotypes this is the same
    # pairwise-diversity quantity underlying Pixy's pi.
    site_pi <- (
        2 * k * (n_chr - k)
    ) / (
        n_chr * (n_chr - 1)
    )

    # --------------------------------------------------------
    # Whole chromosome estimate
    # Kept as a diagnostic only.
    # --------------------------------------------------------

    pi_whole_chr <- sum(site_pi) / CHROM_LEN

    # --------------------------------------------------------
    # 50 kb non-overlapping windows
    # --------------------------------------------------------

    window_numerator <- numeric(n_windows)

    sums <- rowsum(
        site_pi,
        group = variant_window,
        reorder = FALSE
    )

    present_windows <- as.integer(rownames(sums))

    window_numerator[present_windows] <- sums[, 1]

    # All simulated positions are callable.
    # Positions absent from VCF contribute 0 differences but
    # remain part of the window denominator.
    window_pi <- window_numerator / window_length

    # Match Nolen et al.: genome/chromosome mean across windows.
    pi_50kb_mean <- mean(window_pi)

    segregating_sites <- sum(
        k > 0 & k < n_chr
    )

    res[[j]] <- data.frame(
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
        n_chromosomes = n_chr,
        segregating_sites = segregating_sites,
        n_windows = n_windows,
        window_size = WINDOW_SIZE,
        pi_50kb_mean = pi_50kb_mean,
        pi_whole_chr = pi_whole_chr,
        stringsAsFactors = FALSE
    )

    j <- j + 1L
}

out <- do.call(rbind, res)

out <- out[
    order(out$target, out$draw),
    ,
    drop = FALSE
]

# ============================================================
# WRITE
# ============================================================

dir.create(
    "results/pi_pixylike/per_snapshot",
    recursive = TRUE,
    showWarnings = FALSE
)

outfile <- file.path(
    "results/pi_pixylike/per_snapshot",
    sprintf("pi_task%03d.tsv", task_id)
)

write.table(
    out,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Draws calculated       :", nrow(out), "\n")
cat("Mean 50kb-window pi    :", mean(out$pi_50kb_mean), "\n")
cat("Mean whole-chrom pi    :", mean(out$pi_whole_chr), "\n")
cat("Min 50kb-window pi     :", min(out$pi_50kb_mean), "\n")
cat("Max 50kb-window pi     :", max(out$pi_50kb_mean), "\n")
cat("Written                :", outfile, "\n")
