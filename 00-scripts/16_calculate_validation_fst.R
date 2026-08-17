options(stringsAsFactors = FALSE)

CHROM_LEN <- 19361589L
WINDOW    <- 50000L

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1)
    stop("Usage: Rscript 16_calculate_validation_fst.R JOB_ID")

job_id <- as.integer(args[1])

# ============================================================
# 1. READ TASK TABLE + DEFINE 30 INDEPENDENT MODEL REPLICATES
# ============================================================

tasks <- read.delim(
    "results/genotype_export_tasks.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

combos <- unique(
    tasks[, c(
        "landscape",
        "model",
        "rep"
    )]
)

combos <- combos[
    order(
        combos$landscape,
        combos$model,
        combos$rep
    ),
    ,
    drop = FALSE
]

rownames(combos) <- NULL

if (
    job_id < 1 ||
    job_id > nrow(combos)
)
    stop(
        "JOB_ID must be 1-",
        nrow(combos)
    )

combo <- combos[job_id, , drop = FALSE]

LANDSCAPE <- combo$landscape
MODEL     <- combo$model
REP       <- combo$rep

cat("============================================\n")
cat("VALIDATION HUDSON FST\n")
cat("============================================\n")
cat("Job       :", job_id, "\n")
cat("Landscape :", LANDSCAPE, "\n")
cat("Model     :", MODEL, "\n")
cat("Rep       :", REP, "\n")
cat("Window    :", WINDOW, "bp\n\n")

# ============================================================
# 2. COMPARISONS
#
# Historical:
#   E Skane 1956 vs W Skane 1951
#
# Modern:
#   E vs W
#   E vs SE
#   W vs SE
#
# These correspond directly to Nolen et al. Table S7.
# ============================================================

comparisons <- data.frame(

    comparison = c(
        "historical_EW",
        "modern_EW",
        "modern_E_SE",
        "modern_W_SE"
    ),

    target1 = c(
        "ESkane_hist",
        "ESkane_mod",
        "ESkane_mod",
        "WSkane_mod"
    ),

    year1 = c(
        1956,
        2020,
        2020,
        2020
    ),

    target2 = c(
        "WSkane_hist",
        "WSkane_mod",
        "SESkane_mod",
        "SESkane_mod"
    ),

    year2 = c(
        1951,
        2020,
        2020,
        2020
    ),

    empirical_FST = c(
        0.02874,
        0.04940,
        0.08680,
        0.10919
    ),

    stringsAsFactors = FALSE
)

# ============================================================
# 3. READ SAMPLING MANIFEST
# ============================================================

manifest <- read.delim(
    "results/validation_sampling_manifest.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

manifest <- manifest[
    manifest$landscape == LANDSCAPE &
    manifest$model == MODEL &
    manifest$rep == REP,
    ,
    drop = FALSE
]

# ============================================================
# 4. VCF READER
# ============================================================

get_vcf_path <- function(year) {

    x <- tasks[
        tasks$landscape == LANDSCAPE &
        tasks$model == MODEL &
        tasks$rep == REP &
        tasks$year == year,
        ,
        drop = FALSE
    ]

    if (nrow(x) != 1)
        stop(
            "Could not uniquely identify task for ",
            LANDSCAPE, " ",
            MODEL, " ",
            REP, " year ", year,
            ". Found ", nrow(x)
        )

    file.path(
        "validation_genotypes",
        LANDSCAPE,
        MODEL,
        REP,
        paste0(
            LANDSCAPE, "_",
            MODEL, "_",
            REP,
            "_year",
            year,
            ".vcf"
        )
    )
}

read_validation_vcf <- function(path) {

    if (!file.exists(path))
        stop("VCF not found: ", path)

    cat("Reading:", path, "\n")

    con <- file(path, open = "r")

    n_before_header <- 0L

    repeat {

        ln <- readLines(
            con,
            n = 1,
            warn = FALSE
        )

        if (length(ln) == 0)
            stop("No #CHROM header found in ", path)

        if (startsWith(ln, "#CHROM"))
            break

        n_before_header <-
            n_before_header + 1L
    }

    close(con)

    v <- read.delim(
        path,
        skip = n_before_header,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (!"POS" %in% names(v))
        stop("POS column missing in ", path)

    if (ncol(v) < 10)
        stop("No genotype columns in ", path)

    pos <- as.integer(v$POS)

    if (anyDuplicated(pos))
        stop("Duplicate positions found in ", path)

    samples <- names(v)[10:ncol(v)]

    gt <- as.matrix(
        v[, samples, drop = FALSE]
    )

    # Keep only GT if FORMAT ever contains extra fields.
    flat <- sub(
        ":.*$",
        "",
        as.vector(gt)
    )

    alt <- rep(
        NA_integer_,
        length(flat)
    )

    alt[flat %in% c("0|0", "0/0")] <- 0L

    alt[flat %in% c(
        "0|1",
        "1|0",
        "0/1",
        "1/0"
    )] <- 1L

    alt[flat %in% c("1|1", "1/1")] <- 2L

    if (anyNA(alt)) {

        bad <- unique(
            flat[is.na(alt)]
        )

        stop(
            "Unexpected genotype(s): ",
            paste(
                head(bad, 20),
                collapse = ", "
            )
        )
    }

    geno <- matrix(
        alt,
        nrow = nrow(gt),
        ncol = ncol(gt),
        dimnames = list(
            NULL,
            samples
        )
    )

    cat(
        "  variants:", nrow(geno),
        " samples:", ncol(geno),
        "\n"
    )

    list(
        path = path,
        pos = pos,
        geno = geno
    )
}

# ============================================================
# 5. LOAD ONLY REQUIRED SNAPSHOTS
# ============================================================

required_years <- sort(
    unique(
        c(
            comparisons$year1,
            comparisons$year2
        )
    )
)

vcfs <- list()

for (yr in required_years) {

    vcfs[[as.character(yr)]] <-
        read_validation_vcf(
            get_vcf_path(yr)
        )
}

# ============================================================
# 6. HUDSON FST FOR ONE PAIR OF SAMPLES
#
# For each SNP:
#
# pi1, pi2 = unbiased within-population diversity
# dxy      = between-population divergence
#
# Hudson numerator:
#   dxy - (pi1 + pi2)/2
#
# Hudson denominator:
#   dxy
#
# Components are summed within 50-kb windows.
# Mean FST is then arithmetic mean across windows,
# matching the empirical analysis.
# ============================================================

calculate_fst <- function(
    v1,
    inds1,
    v2,
    inds2
) {

    s1 <- paste0(
        "p1:i",
        inds1
    )

    s2 <- paste0(
        "p1:i",
        inds2
    )

    idx1 <- match(
        s1,
        colnames(v1$geno)
    )

    idx2 <- match(
        s2,
        colnames(v2$geno)
    )

    if (anyNA(idx1))
        stop(
            "Population 1 sample missing from VCF."
        )

    if (anyNA(idx2))
        stop(
            "Population 2 sample missing from VCF."
        )

    k1_raw <- rowSums(
        v1$geno[
            ,
            idx1,
            drop = FALSE
        ]
    )

    k2_raw <- rowSums(
        v2$geno[
            ,
            idx2,
            drop = FALSE
        ]
    )

    # --------------------------------------------------------
    # Align positions.
    #
    # For comparisons involving different years, use the union
    # of variant positions. If a mutation is absent from one
    # snapshot, all sampled chromosomes there carry reference.
    # --------------------------------------------------------

    if (
        identical(
            v1$path,
            v2$path
        )
    ) {

        pos <- v1$pos
        k1  <- k1_raw
        k2  <- k2_raw

    } else {

        pos <- sort(
            unique(
                c(
                    v1$pos,
                    v2$pos
                )
            )
        )

        k1 <- integer(
            length(pos)
        )

        k2 <- integer(
            length(pos)
        )

        k1[
            match(
                v1$pos,
                pos
            )
        ] <- k1_raw

        k2[
            match(
                v2$pos,
                pos
            )
        ] <- k2_raw
    }

    n1 <- 2L * length(inds1)
    n2 <- 2L * length(inds2)

    # SNPs segregating across the two sampled populations.
    variable <- (
        (k1 + k2) > 0 &
        (k1 + k2) < (n1 + n2)
    )

    pos <- pos[variable]
    k1  <- k1[variable]
    k2  <- k2[variable]

    if (length(pos) == 0)
        return(
            c(
                mean_window_FST = NA,
                component_FST = NA,
                n_windows = 0,
                n_snps = 0
            )
        )

    # Unbiased within-population nucleotide diversity.
    pi1 <- (
        2 * k1 * (n1 - k1) /
        (n1 * (n1 - 1))
    )

    pi2 <- (
        2 * k2 * (n2 - k2) /
        (n2 * (n2 - 1))
    )

    # Between-population pairwise divergence.
    dxy <- (
        k1 * (n2 - k2) +
        (n1 - k1) * k2
    ) / (
        n1 * n2
    )

    numerator <- (
        dxy -
        (pi1 + pi2) / 2
    )

    denominator <- dxy

    window <- (
        (pos - 1L) %/% WINDOW
    ) + 1L

    win_num <- tapply(
        numerator,
        window,
        sum
    )

    win_den <- tapply(
        denominator,
        window,
        sum
    )

    win_fst <- (
        win_num /
        win_den
    )

    valid <- (
        is.finite(win_fst) &
        win_den > 0
    )

    # Primary statistic:
    # arithmetic mean across 50-kb window estimates,
    # matching the approach described by Nolen et al.
    mean_window_FST <- mean(
        win_fst[valid]
    )

    # Additional QC statistic:
    # ratio of total Hudson components across chromosome.
    component_FST <- (
        sum(numerator) /
        sum(denominator)
    )

    c(
        mean_window_FST =
            mean_window_FST,

        component_FST =
            component_FST,

        n_windows =
            sum(valid),

        n_snps =
            length(pos)
    )
}

# ============================================================
# 7. 100 EMPIRICAL-SIZED DRAWS PER COMPARISON
# ============================================================

res <- list()
z <- 1L

for (cc in seq_len(nrow(comparisons))) {

    cmp <- comparisons[cc, ]

    cat(
        "\nComparison:",
        cmp$comparison,
        "\n"
    )

    m1 <- manifest[
        manifest$target == cmp$target1 &
        manifest$simulation_year == cmp$year1,
        ,
        drop = FALSE
    ]

    m2 <- manifest[
        manifest$target == cmp$target2 &
        manifest$simulation_year == cmp$year2,
        ,
        drop = FALSE
    ]

    if (
        nrow(m1) == 0 ||
        nrow(m2) == 0
    )
        stop(
            "Missing manifest data for ",
            cmp$comparison
        )

    draws <- intersect(
        unique(m1$draw),
        unique(m2$draw)
    )

    draws <- sort(draws)

    if (length(draws) != 100)
        stop(
            cmp$comparison,
            ": expected 100 shared draws, found ",
            length(draws)
        )

    v1 <- vcfs[[as.character(cmp$year1)]]

    v2 <- vcfs[[as.character(cmp$year2)]]

    for (dd in draws) {

        g1 <- m1[
            m1$draw == dd,
            ,
            drop = FALSE
        ]

        g2 <- m2[
            m2$draw == dd,
            ,
            drop = FALSE
        ]

        if (anyDuplicated(g1$ind_index))
            stop(
                "Duplicate population-1 individual in ",
                cmp$comparison, " ", dd
            )

        if (anyDuplicated(g2$ind_index))
            stop(
                "Duplicate population-2 individual in ",
                cmp$comparison, " ", dd
            )

        fst <- calculate_fst(
            v1,
            g1$ind_index,
            v2,
            g2$ind_index
        )

        res[[z]] <- data.frame(

            job_id = job_id,

            landscape = LANDSCAPE,
            model = MODEL,
            rep = REP,

            comparison =
                cmp$comparison,

            target1 =
                cmp$target1,

            year1 =
                cmp$year1,

            n1 =
                nrow(g1),

            target2 =
                cmp$target2,

            year2 =
                cmp$year2,

            n2 =
                nrow(g2),

            draw = dd,

            mean_window_FST =
                fst[
                    "mean_window_FST"
                ],

            component_FST =
                fst[
                    "component_FST"
                ],

            n_windows =
                as.integer(
                    fst["n_windows"]
                ),

            n_snps =
                as.integer(
                    fst["n_snps"]
                ),

            empirical_FST =
                cmp$empirical_FST,

            error =
                fst[
                    "mean_window_FST"
                ] -
                cmp$empirical_FST,

            stringsAsFactors = FALSE
        )

        z <- z + 1L
    }
}

out <- do.call(
    rbind,
    res
)

rownames(out) <- NULL

# ============================================================
# 8. WRITE DRAW-LEVEL OUTPUT
# ============================================================

outfile <- file.path(
    "results/fst/per_modelrep",
    sprintf(
        "FST_job%03d_draws.tsv",
        job_id
    )
)

write.table(
    out,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 9. QC SUMMARY
# ============================================================

cat("\n============================================\n")
cat("FST QC SUMMARY\n")
cat("============================================\n")

for (
    cmp_name in
    unique(out$comparison)
) {

    x <- out[
        out$comparison ==
            cmp_name,
        ,
        drop = FALSE
    ]

    cat("\n", cmp_name, "\n", sep = "")

    cat(
        "  empirical FST :",
        unique(x$empirical_FST),
        "\n"
    )

    cat(
        "  simulation mean:",
        mean(x$mean_window_FST),
        "\n"
    )

    cat(
        "  simulation SD  :",
        sd(x$mean_window_FST),
        "\n"
    )

    cat(
        "  range          :",
        min(x$mean_window_FST),
        "-",
        max(x$mean_window_FST),
        "\n"
    )

    cat(
        "  mean windows   :",
        mean(x$n_windows),
        "\n"
    )

    cat(
        "  mean SNPs      :",
        mean(x$n_snps),
        "\n"
    )

    cat(
        "  component FST  :",
        mean(x$component_FST),
        "\n"
    )
}

# Direct historical -> modern E-W increase
hist <- out[
    out$comparison ==
        "historical_EW",
]

modern <- out[
    out$comparison ==
        "modern_EW",
]

sim_delta <- (
    mean(modern$mean_window_FST) -
    mean(hist$mean_window_FST)
)

emp_delta <- (
    0.04940 -
    0.02874
)

cat("\nE-W temporal change:\n")
cat(
    "  simulated delta:",
    sim_delta,
    "\n"
)
cat(
    "  empirical delta:",
    emp_delta,
    "\n"
)
cat(
    "  error          :",
    sim_delta - emp_delta,
    "\n"
)

cat("\nWritten:\n")
cat(outfile, "\n")
