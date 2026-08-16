# ============================================================
# Swedish Butterfly Monitoring Scheme density summary
#
# Inputs:
#   01-data/monitoring/semiargus_densities_along_transects.tsv
#   01-data/monitoring/semiargus_densities_along_transect_segments.tsv
#   01-data/monitoring/semiargus_densities_at_point_sites.tsv
#   01-data/monitoring/corine_index_to_clc.tsv
#   02-landscapes/derived/segments_annotated.tsv
#   02-landscapes/derived/point_densities_with_corine.tsv
#
# Outputs:
#   02-landscapes/derived/density_records_standardized.tsv
#   02-landscapes/tables/density_summary_by_method.tsv
#   02-landscapes/tables/density_summary_by_habitat.tsv
#   02-landscapes/figures/density_by_habitat.*
#
# Usage:
#   Rscript 00-scripts/03d_density_summary.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ============================================================
# 1. PATHS
# ============================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Run this script using Rscript.")
}

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))

MON <- file.path(REPO_ROOT, "01-data", "monitoring")
DER <- file.path(REPO_ROOT, "02-landscapes", "derived")
TAB <- file.path(REPO_ROOT, "02-landscapes", "tables")
FIG <- file.path(REPO_ROOT, "02-landscapes", "figures")

dir.create(DER, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# Pollard walks count within 2.5 m on either side of the route.
TRANSECT_WIDTH_M <- 5

# 1000 m × 5 m = 5000 m² = 0.005 km².
# Therefore counts per 1000 m × 200 = observed counts per km².
PER_1000M_TO_KM2 <- 1e6 / (1000 * TRANSECT_WIDTH_M)

# ============================================================
# 2. INPUT DATA
# ============================================================

transects <- read_tsv(
  file.path(MON, "semiargus_densities_along_transects.tsv"),
  show_col_types = FALSE
)

segments <- read_tsv(
  file.path(MON, "semiargus_densities_along_transect_segments.tsv"),
  show_col_types = FALSE
)

points <- read_tsv(
  file.path(MON, "semiargus_densities_at_point_sites.tsv"),
  show_col_types = FALSE
)

segments_ann <- read_tsv(
  file.path(DER, "segments_annotated.tsv"),
  show_col_types = FALSE
)

points_ann <- read_tsv(
  file.path(DER, "point_densities_with_corine.tsv"),
  show_col_types = FALSE
)

clc_lookup <- read_tsv(
  file.path(MON, "corine_index_to_clc.tsv"),
  show_col_types = FALSE
) |>
  select(corine1, clc_code, clc_name) |>
  distinct()

# ============================================================
# 3. CHECK THE ORIGINAL STANDARDISATION
# ============================================================

cat("\n==================== STANDARDISATION CHECKS ====================\n")

# Point sites cover a circle of radius 25 m.
point_area_m2 <- pi * 25^2
expected_point_equivalent_length <- point_area_m2 / TRANSECT_WIDTH_M

cat(
  "Point-site survey area:", round(point_area_m2, 1), "m2\n"
)

cat(
  "Equivalent 5-m transect length:",
  round(expected_point_equivalent_length, 1), "m\n"
)

cat(
  "Median sit_length_equivalent in data:",
  round(median(points$sit_length_equivalent, na.rm = TRUE), 1),
  "m\n"
)

# Check standardized count calculations.
transect_error <- max(
  abs(
    transects$count_per_1000m -
      transects$sum / transects$length_of_transect * 1000
  ),
  na.rm = TRUE
)

segment_error <- max(
  abs(
    segments$count_per_1000m -
      segments$sum / segments$segment_length * 1000
  ),
  na.rm = TRUE
)

point_error <- max(
  abs(
    points$count_per_1000m_equivalent -
      points$sum / points$sit_length_equivalent * 1000
  ),
  na.rm = TRUE
)

cat("Maximum transect scaling error:", signif(transect_error, 3), "\n")
cat("Maximum segment scaling error:", signif(segment_error, 3), "\n")
cat("Maximum point scaling error:", signif(point_error, 3), "\n")

# ============================================================
# 4. STANDARDIZE TO OBSERVED COUNTS PER KM2
# ============================================================

transect_std <- transects |>
  transmute(
    method = "Whole transect",
    sit_uid,
    sit_name,
    year,
    count = sum,
    count_per_1000m,
    observed_density_km2 =
      count_per_1000m * PER_1000M_TO_KM2
  )

segment_std <- segments |>
  transmute(
    method = "Transect segment",
    sit_uid,
    sit_name,
    year,
    count = sum,
    count_per_1000m,
    observed_density_km2 =
      count_per_1000m * PER_1000M_TO_KM2
  )

point_std <- points |>
  transmute(
    method = "Point site",
    sit_uid,
    sit_name,
    year,
    count = sum,
    count_per_1000m = count_per_1000m_equivalent,
    observed_density_km2 =
      count_per_1000m_equivalent * PER_1000M_TO_KM2
  )

density_all <- bind_rows(
  transect_std,
  segment_std,
  point_std
)

write_tsv(
  density_all,
  file.path(DER, "density_records_standardized.tsv")
)

# ============================================================
# 5. OVERALL DENSITY RANGE
# ============================================================

density_summary <- density_all |>
  group_by(method) |>
  summarise(
    n_records = n(),
    n_sites = n_distinct(sit_uid),
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),

    n_zero = sum(count == 0, na.rm = TRUE),

    min = min(observed_density_km2, na.rm = TRUE),
    q05 = quantile(observed_density_km2, 0.05, na.rm = TRUE),
    q25 = quantile(observed_density_km2, 0.25, na.rm = TRUE),
    median = median(observed_density_km2, na.rm = TRUE),
    mean = mean(observed_density_km2, na.rm = TRUE),
    q75 = quantile(observed_density_km2, 0.75, na.rm = TRUE),
    q95 = quantile(observed_density_km2, 0.95, na.rm = TRUE),
    max = max(observed_density_km2, na.rm = TRUE),

    .groups = "drop"
  )

write_tsv(
  density_summary,
  file.path(TAB, "density_summary_by_method.tsv")
)

cat("\n==================== OBSERVED DENSITY RANGE ====================\n")

print(
  as.data.frame(
    density_summary |>
      mutate(
        across(
          c(min, q05, q25, median, mean, q75, q95, max),
          \(x) round(x, 1)
        )
      )
  )
)

if (all(density_summary$n_zero == 0)) {
  cat(
    "\nNOTE: No zero-count observations occur in these density files.\n",
    "The summaries therefore describe records with observed C. semiargus,\n",
    "not average density across all SeBMS survey visits.\n",
    sep = ""
  )
}

# ============================================================
# 6. HABITAT CLASSIFICATION
# ============================================================

clc_to_group <- function(code) {

  code <- as.integer(code)

  case_when(
    is.na(code) ~ NA_character_,
    code == 231 ~ "Grassland",
    code %in% c(321, 322, 323, 324) ~ "Grassland",

    code %in% c(
      211, 212, 213,
      221, 222, 223,
      241, 242, 243, 244
    ) ~ "Cropland",

    code %in% c(311, 312, 313) ~ "Forest",

    code >= 111 & code <= 142 ~ "Urban",

    TRUE ~ "Other"
  )
}

# QGIS annotation stores the CORINE raster category index.
# Convert that index to the actual CLC code using corine_index_to_clc.tsv.
attach_corine <- function(d) {

  possible_index_cols <- c(
    "corine_1",
    "corine1",
    "clc_index",
    "clc_code"
  )

  found <- intersect(possible_index_cols, names(d))

  if (length(found) == 0) {
    stop(
      "No CORINE index/code column found. Columns are: ",
      paste(names(d), collapse = ", ")
    )
  }

  col <- found[1]

  if (col == "clc_code") {

    d |>
      mutate(
        clc_code = as.integer(clc_code),
        habitat = clc_to_group(clc_code)
      )

  } else {

    d |>
      mutate(
        corine1 = as.integer(.data[[col]])
      ) |>
      left_join(
        clc_lookup,
        by = "corine1"
      ) |>
      mutate(
        habitat = clc_to_group(clc_code)
      )
  }
}

segments_hab <- attach_corine(segments_ann) |>
  transmute(
    method = "Transect segment",
    sit_uid,
    sit_name,
    year,
    count = sum,
    count_per_1000m,
    observed_density_km2 =
      count_per_1000m * PER_1000M_TO_KM2,
    clc_code,
    clc_name,
    habitat
  )

points_hab_raw <- attach_corine(points_ann)

# Point file may retain either the original name or the standardized name.
if ("count_per_1000m_equivalent" %in% names(points_hab_raw)) {

  points_hab <- points_hab_raw |>
    transmute(
      method = "Point site",
      sit_uid,
      sit_name,
      year,
      count = sum,
      count_per_1000m = count_per_1000m_equivalent,
      observed_density_km2 =
        count_per_1000m_equivalent * PER_1000M_TO_KM2,
      clc_code,
      clc_name,
      habitat
    )

} else if ("count_per_1000m" %in% names(points_hab_raw)) {

  points_hab <- points_hab_raw |>
    transmute(
      method = "Point site",
      sit_uid,
      sit_name,
      year,
      count = sum,
      count_per_1000m,
      observed_density_km2 =
        count_per_1000m * PER_1000M_TO_KM2,
      clc_code,
      clc_name,
      habitat
    )

} else {

  stop(
    "Point annotation file contains neither ",
    "count_per_1000m_equivalent nor count_per_1000m."
  )
}

habitat_density <- bind_rows(
  segments_hab,
  points_hab
) |>
  filter(!is.na(habitat))

# ============================================================
# 7. DENSITY BY HABITAT
# ============================================================

habitat_summary <- habitat_density |>
  group_by(method, habitat) |>
  summarise(
    n_records = n(),
    n_sites = n_distinct(sit_uid),

    q05 = quantile(observed_density_km2, 0.05, na.rm = TRUE),
    q25 = quantile(observed_density_km2, 0.25, na.rm = TRUE),
    median = median(observed_density_km2, na.rm = TRUE),
    mean = mean(observed_density_km2, na.rm = TRUE),
    q75 = quantile(observed_density_km2, 0.75, na.rm = TRUE),
    q95 = quantile(observed_density_km2, 0.95, na.rm = TRUE),
    max = max(observed_density_km2, na.rm = TRUE),

    .groups = "drop"
  )

write_tsv(
  habitat_summary,
  file.path(TAB, "density_summary_by_habitat.tsv")
)

cat("\n==================== DENSITY BY HABITAT ====================\n")

print(
  as.data.frame(
    habitat_summary |>
      mutate(
        across(
          c(q05, q25, median, mean, q75, q95, max),
          \(x) round(x, 1)
        )
      )
  )
)

# ============================================================
# 8. FIGURE
# ============================================================

HABITAT_ORDER <- c(
  "Grassland",
  "Cropland",
  "Forest",
  "Urban",
  "Other"
)

plot_dat <- habitat_density |>
  mutate(
    habitat = factor(habitat, levels = HABITAT_ORDER)
  )

p <- ggplot(
  plot_dat,
  aes(
    x = habitat,
    y = observed_density_km2
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 1
  ) +
  facet_wrap(
    ~method,
    scales = "free_y"
  ) +
  scale_y_log10() +
  labs(
    x = "CORINE habitat class",
    y = expression(
      "Area-standardized observed density (butterflies " * km^{-2} * ")"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  file.path(FIG, "density_by_habitat.png"),
  p,
  width = 8,
  height = 4.5,
  dpi = 320
)

ggsave(
  file.path(FIG, "density_by_habitat.pdf"),
  p,
  width = 8,
  height = 4.5
)

ggsave(
  file.path(FIG, "density_by_habitat.svg"),
  p,
  width = 8,
  height = 4.5
)

# ============================================================
# 9. OUTPUT SUMMARY
# ============================================================

cat("\n==================== OUTPUTS ====================\n")

cat(
  "Overall density range:\n",
  "  density_summary_by_method.tsv\n\n",
  sep = ""
)

cat(
  "Density by CORINE habitat:\n",
  "  density_summary_by_habitat.tsv\n",
  "  density_by_habitat.*\n\n",
  sep = ""
)

cat(
  "Standardized records:\n",
  "  density_records_standardized.tsv\n",
  sep = ""
)
