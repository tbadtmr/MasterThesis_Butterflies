options(stringsAsFactors = FALSE)

CHROM_LEN <- 19361589
MIN_ROH   <- 100000

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        "Usage: Rscript 26_parse_regional_roh.R ",
        "RAW_ROH MANIFEST OUTROOT BASE"
    )
}

raw_file     <- args[1]
manifest_file <- args[2]
outroot      <- args[3]
base         <- args[4]

# ============================================================
# READ MANIFEST
# ============================================================

manifest <- read.delim(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!"rep" %in% names(manifest) &&
    "replicate" %in% names(manifest)) {
    manifest$rep <- manifest$replicate
}

required <- c(
    "scenario",
    "rep",
    "year",
    "region",
    "site",
    "site_number",
    "ind_index"
)

missing_cols <- setdiff(
    required,
    names(manifest)
)

if (length(missing_cols) > 0) {
    stop(
        "Manifest missing columns: ",
        paste(missing_cols, collapse = ", ")
    )
}

if (nrow(manifest) == 0) {
    stop("Filtered site-2 manifest is empty.")
}

if (anyDuplicated(manifest$ind_index)) {
    stop("Duplicate ind_index in regional manifest.")
}

scenario <- unique(manifest$scenario)
rep_id   <- unique(manifest$rep)
year     <- unique(manifest$year)

if (
    length(scenario) != 1 ||
    length(rep_id) != 1 ||
    length(year) != 1
) {
    stop("Manifest does not contain one unique scenario/rep/year.")
}

# Site-2 regional sampling QC
region_counts <- table(manifest$region)

if (
    length(region_counts) != 5 ||
    any(region_counts < 1) ||
    any(region_counts > 20)
) {
    stop(
        "Expected five regions with 1-20 individuals each. Found:
",
        paste(
            names(region_counts),
            region_counts,
            collapse = "
"
        )
    )
}

if (any(manifest$site_number != 2)) {
    stop("Manifest contains individuals outside site_number = 2.")
}

sites_per_region <- tapply(
    manifest$site,
    manifest$region,
    function(z) length(unique(z))
)

if (any(sites_per_region != 1)) {
    stop("More than one sampling site found within a region.")
}

manifest$sample_name <- paste0(
    "p1:i",
    manifest$ind_index
)

# ============================================================
# READ RAW BCFTOOLS/ROH OUTPUT
# ============================================================

lines <- readLines(
    raw_file,
    warn = FALSE
)

rg_lines <- lines[
    grepl("^RG\\t", lines)
]

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
        text = paste(
            rg_lines,
            collapse = "\n"
        ),
        header = FALSE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE
    )

    if (ncol(roh) != 8) {
        stop(
            "Unexpected number of ROH columns: ",
            ncol(roh)
        )
    }

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
cat("REGIONAL ROH ANALYSIS\n")
cat("============================================\n")
cat("Scenario        :", scenario, "\n")
cat("Replicate       :", rep_id, "\n")
cat("Year            :", year, "\n")
cat("Individuals     :", nrow(manifest), "\n")
cat("Chromosome len  :", CHROM_LEN, "\n")
cat("All RG segments :", nrow(roh), "\n")

# ============================================================
# RETAIN ROH >100 kb
# ============================================================

roh100 <- roh[
    roh$length_bp > MIN_ROH,
    ,
    drop = FALSE
]

cat("ROH >100 kb     :", nrow(roh100), "\n")

# ============================================================
# CHECK SAMPLE NAMES
# ============================================================

if (nrow(roh100) > 0) {

    unknown <- setdiff(
        unique(roh100$sample_name),
        manifest$sample_name
    )

    if (length(unknown) > 0) {
        stop(
            "ROH output contains unknown sample names: ",
            paste(
                head(unknown, 20),
                collapse = ", "
            )
        )
    }
}

# ============================================================
# ROH LENGTH CLASSES
#
# Keep raw length too, so bins can always be changed later.
# ============================================================

if (nrow(roh100) > 0) {

    roh100$roh_class <- cut(
        roh100$length_bp,
        breaks = c(
            MIN_ROH,
            250000,
            500000,
            1000000,
            Inf
        ),
        labels = c(
            "100-250kb",
            "250-500kb",
            "500kb-1Mb",
            ">1Mb"
        ),
        right = TRUE,
        include.lowest = FALSE
    )

    ii <- match(
        roh100$sample_name,
        manifest$sample_name
    )

    roh100$scenario <- scenario
    roh100$rep <- rep_id
    roh100$year <- year

    roh100$region <- manifest$region[ii]
    roh100$site <- manifest$site[ii]
    roh100$ind_index <- manifest$ind_index[ii]

    # Put metadata first
    roh100 <- roh100[, c(
        "scenario",
        "rep",
        "year",
        "region",
        "site",
        "sample_name",
        "ind_index",
        "chromosome",
        "start",
        "end",
        "length_bp",
        "roh_class",
        "n_markers",
        "quality"
    )]

} else {

    roh100 <- data.frame(
        scenario = character(),
        rep = character(),
        year = integer(),
        region = character(),
        site = character(),
        sample_name = character(),
        ind_index = integer(),
        chromosome = character(),
        start = integer(),
        end = integer(),
        length_bp = integer(),
        roh_class = character(),
        n_markers = integer(),
        quality = numeric(),
        stringsAsFactors = FALSE
    )
}

# ============================================================
# INDIVIDUAL-LEVEL TABLE
# ============================================================

individual <- manifest[, c(
    "scenario",
    "rep",
    "year",
    "region",
    "site",
    "site_x",
    "site_y",
    "sample_name",
    "ind_index"
)]

individual$n_roh_gt100kb <- 0L
individual$total_roh_bp <- 0
individual$mean_roh_bp <- NA_real_
individual$median_roh_bp <- NA_real_
individual$max_roh_bp <- 0

individual$n_roh_100_250kb <- 0L
individual$n_roh_250_500kb <- 0L
individual$n_roh_500kb_1Mb <- 0L
individual$n_roh_gt1Mb <- 0L

individual$bp_roh_100_250kb <- 0
individual$bp_roh_250_500kb <- 0
individual$bp_roh_500kb_1Mb <- 0
individual$bp_roh_gt1Mb <- 0

if (nrow(roh100) > 0) {

    groups <- split(
        roh100,
        roh100$sample_name
    )

    for (s in names(groups)) {

        g <- groups[[s]]

        ii <- match(
            s,
            individual$sample_name
        )

        if (is.na(ii)) {
            stop(
                "Could not match sample ",
                s
            )
        }

        individual$n_roh_gt100kb[ii] <- nrow(g)

        individual$total_roh_bp[ii] <- sum(
            g$length_bp
        )

        individual$mean_roh_bp[ii] <- mean(
            g$length_bp
        )

        individual$median_roh_bp[ii] <- median(
            g$length_bp
        )

        individual$max_roh_bp[ii] <- max(
            g$length_bp
        )

        # --------------------------------------------
        # Length-class contributions
        # --------------------------------------------

        cls <- as.character(
            g$roh_class
        )

        individual$n_roh_100_250kb[ii] <- sum(
            cls == "100-250kb"
        )

        individual$n_roh_250_500kb[ii] <- sum(
            cls == "250-500kb"
        )

        individual$n_roh_500kb_1Mb[ii] <- sum(
            cls == "500kb-1Mb"
        )

        individual$n_roh_gt1Mb[ii] <- sum(
            cls == ">1Mb"
        )

        individual$bp_roh_100_250kb[ii] <- sum(
            g$length_bp[
                cls == "100-250kb"
            ]
        )

        individual$bp_roh_250_500kb[ii] <- sum(
            g$length_bp[
                cls == "250-500kb"
            ]
        )

        individual$bp_roh_500kb_1Mb[ii] <- sum(
            g$length_bp[
                cls == "500kb-1Mb"
            ]
        )

        individual$bp_roh_gt1Mb[ii] <- sum(
            g$length_bp[
                cls == ">1Mb"
            ]
        )
    }
}

# ============================================================
# FROH
# ============================================================

individual$FROH_100kb <- (
    individual$total_roh_bp /
    CHROM_LEN
)

individual$FROH_100_250kb <- (
    individual$bp_roh_100_250kb /
    CHROM_LEN
)

individual$FROH_250_500kb <- (
    individual$bp_roh_250_500kb /
    CHROM_LEN
)

individual$FROH_500kb_1Mb <- (
    individual$bp_roh_500kb_1Mb /
    CHROM_LEN
)

individual$FROH_gt1Mb <- (
    individual$bp_roh_gt1Mb /
    CHROM_LEN
)

if (
    any(
        individual$FROH_100kb >
            1 + 1e-8
    )
) {
    stop("FROH > 1 detected.")
}

# Check class decomposition
class_sum <- (
    individual$bp_roh_100_250kb +
    individual$bp_roh_250_500kb +
    individual$bp_roh_500kb_1Mb +
    individual$bp_roh_gt1Mb
)

if (
    any(
        abs(
            class_sum -
            individual$total_roh_bp
        ) > 1
    )
) {
    stop(
        "ROH length-class decomposition failed."
    )
}

# ============================================================
# REGIONAL SUMMARY
#
# IMPORTANT:
# Summaries are first calculated per individual, then averaged
# within region. We do not pool ROH segments directly.
# ============================================================

mean_or_na <- function(x) {

    x <- x[
        is.finite(x)
    ]

    if (length(x) == 0) {
        return(NA_real_)
    }

    mean(x)
}

sd_or_na <- function(x) {

    x <- x[
        is.finite(x)
    ]

    if (length(x) < 2) {
        return(NA_real_)
    }

    sd(x)
}

region_order <- c(
    "NW",
    "NE",
    "W",
    "E",
    "SE"
)

regional_list <- lapply(
    region_order,
    function(reg) {

        g <- individual[
            individual$region == reg,
            ,
            drop = FALSE
        ]

        seg <- roh100[
            roh100$region == reg,
            ,
            drop = FALSE
        ]

        data.frame(
            scenario = scenario,
            rep = rep_id,
            year = year,
            region = reg,

            n_individuals = nrow(g),
            n_individuals_with_roh =
                sum(
                    g$n_roh_gt100kb > 0
                ),

            mean_FROH_100kb =
                mean(
                    g$FROH_100kb
                ),

            sd_FROH_100kb =
                sd_or_na(
                    g$FROH_100kb
                ),

            median_FROH_100kb =
                median(
                    g$FROH_100kb
                ),

            mean_n_roh =
                mean(
                    g$n_roh_gt100kb
                ),

            mean_roh_length_bp =
                mean_or_na(
                    g$mean_roh_bp
                ),

            mean_individual_max_roh_bp =
                mean(
                    g$max_roh_bp
                ),

            max_roh_bp_observed =
                if (nrow(seg) > 0)
                    max(seg$length_bp)
                else
                    0,

            mean_FROH_100_250kb =
                mean(
                    g$FROH_100_250kb
                ),

            mean_FROH_250_500kb =
                mean(
                    g$FROH_250_500kb
                ),

            mean_FROH_500kb_1Mb =
                mean(
                    g$FROH_500kb_1Mb
                ),

            mean_FROH_gt1Mb =
                mean(
                    g$FROH_gt1Mb
                ),

            stringsAsFactors = FALSE
        )
    }
)

regional <- do.call(
    rbind,
    regional_list
)

# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

seg_dir <- file.path(
    outroot,
    "segments",
    scenario
)

ind_dir <- file.path(
    outroot,
    "individuals",
    scenario
)

sum_dir <- file.path(
    outroot,
    "summaries",
    scenario
)

dir.create(
    seg_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    ind_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    sum_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# WRITE OUTPUT
# ============================================================

segment_file <- file.path(
    seg_dir,
    paste0(
        base,
        "_ROH_segments_gt100kb.tsv"
    )
)

individual_file <- file.path(
    ind_dir,
    paste0(
        base,
        "_individual_ROH_metrics.tsv"
    )
)

regional_file <- file.path(
    sum_dir,
    paste0(
        base,
        "_regional_ROH_summary.tsv"
    )
)

write.table(
    roh100,
    segment_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    individual,
    individual_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    regional,
    regional_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# FINAL QC
# ============================================================

if (nrow(individual) != nrow(manifest)) {
    stop(
        "Individual output size does not match manifest."
    )
}

if (nrow(regional) != 5) {
    stop(
        "Regional summary does not contain five rows."
    )
}

cat("\n")
cat("Mean FROH by region:\n")

print(
    regional[, c(
        "region",
        "mean_FROH_100kb",
        "mean_n_roh",
        "max_roh_bp_observed"
    )],
    row.names = FALSE
)

cat("\n")
cat("Segments:   ", segment_file, "\n")
cat("Individuals:", individual_file, "\n")
cat("Regional:   ", regional_file, "\n")

cat("\nDONE\n")
