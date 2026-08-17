options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2)
    stop("Usage: Rscript 15_compare_froh_rates.R DIR_A DIR_B")

read_change <- function(d) {
    f <- file.path("results", d, "FROH_change_by_model.tsv")
    if (!file.exists(f)) stop("Missing: ", f)
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
}

a <- read_change(args[1])
b <- read_change(args[2])

keys <- c("landscape", "model", "comparison")

cols <- c(keys,
          "mean_historical_FROH",
          "mean_modern_FROH",
          "mean_delta_FROH",
          "sd_delta_FROH")

m <- merge(
    a[, c(cols, "empirical_delta_FROH")],
    b[, cols],
    by = keys,
    suffixes = c("_A", "_B")
)

m$ratio_hist  <- m$mean_historical_FROH_B / m$mean_historical_FROH_A
m$ratio_mod   <- m$mean_modern_FROH_B / m$mean_modern_FROH_A
m$ratio_delta <- m$mean_delta_FROH_B / m$mean_delta_FROH_A

m$error_A <- m$mean_delta_FROH_A - m$empirical_delta_FROH
m$error_B <- m$mean_delta_FROH_B - m$empirical_delta_FROH

m <- m[order(m$landscape, m$model, m$comparison), , drop = FALSE]

outfile <- file.path("results", "FROH_rate_comparison.tsv")

write.table(m, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

cat("============================================\n")
cat("FROH RECOMBINATION-RATE COMPARISON\n")
cat("============================================\n")
cat("A (primary):", args[1], "\n")
cat("B (variant):", args[2], "\n\n")

show <- m[, c("landscape", "model", "comparison",
              "mean_historical_FROH_A", "mean_historical_FROH_B",
              "mean_modern_FROH_A", "mean_modern_FROH_B",
              "mean_delta_FROH_A", "mean_delta_FROH_B",
              "empirical_delta_FROH")]

for (j in seq_along(show))
    if (is.numeric(show[[j]]))
        show[[j]] <- signif(show[[j]], 3)

print(show, row.names = FALSE)

cat("\nMedian ratio B/A, historical FROH:",
    signif(median(m$ratio_hist, na.rm = TRUE), 3), "\n")
cat("Median ratio B/A, modern FROH    :",
    signif(median(m$ratio_mod, na.rm = TRUE), 3), "\n")
cat("Median ratio B/A, delta FROH     :",
    signif(median(m$ratio_delta, na.rm = TRUE), 3), "\n")

cat("\nMean absolute error in delta FROH:\n")
cat("  ", args[1], ":", signif(mean(abs(m$error_A)), 4), "\n")
cat("  ", args[2], ":", signif(mean(abs(m$error_B)), 4), "\n")

cat("\nWritten:", outfile, "\n")
