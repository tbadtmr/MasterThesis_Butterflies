# ============================================================
# GBIF habitat association analysis
#
# Inputs:
#   01-data/gbif/occurrence.txt
#   01-data/monitoring/corine_index_to_clc.tsv
#   01-data/CORINE/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif
#   02-landscapes/derived/sites_with_classes.tsv
#
# Outputs:
#   02-landscapes/tables/TABLE_S4_gbif_habitat_selection.tsv
#   02-landscapes/tables/TABLE_S4_gbif_habitat_selection_numeric.tsv
#   02-landscapes/tables/gbif_record_filtering.tsv
#
# Usage:
#   Rscript 00-scripts/03e_gbif_habitat_analysis.R
#
# GBIF download:
#   GBIF.org (17 March 2026), DOI: 10.15468/dl.9kxv2n
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
})

set.seed(123)

# ============================================================
# 1. PATHS AND SETTINGS
# ============================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Run this script using Rscript.")
}

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))

GBIF_IN <- file.path(
  REPO_ROOT,
  "01-data", "gbif", "occurrence.txt"
)

CLC_LOOKUP <- file.path(
  REPO_ROOT,
  "01-data", "monitoring", "corine_index_to_clc.tsv"
)

SITES_IN <- file.path(
  REPO_ROOT,
  "02-landscapes", "derived", "sites_with_classes.tsv"
)

CORINE_TIF <- file.path(
  REPO_ROOT,
  "01-data", "CORINE",
  "u2018_clc2018_v2020_20u1_raster100m",
  "DATA",
  "U2018_CLC2018_V2020_20u1.tif"
)

TAB <- file.path(
  REPO_ROOT,
  "02-landscapes", "tables"
)

dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

SITE_CRS <- "EPSG:3006"

# Simulation study-region bounds in SWEREF99 TM.
STUDY_YMIN <- 6112797
STUDY_YMAX <- 6418978
STUDY_XMIN <- 295217
STUDY_XMAX <- 672595

# Analysis filters.
MAX_COORD_UNCERTAINTY <- 1000
FLIGHT_MONTHS <- 5:8
YEAR_MIN <- 2000
YEAR_MAX <- 2025

# Spatial thinning and background sampling.
THIN_GRID_M <- 1000
N_BOOT       <- 999

GROUP_ORDER <- c(
  "Grassland",
  "Cropland",
  "Forest",
  "Urban",
  "Other"
)

# ============================================================
# 2. HABITAT GROUPING
# ============================================================

# Same broad CORINE grouping used for the SeBMS analyses.

clc_to_group <- function(code) {

  code <- as.integer(code)

  case_when(
    is.na(code) ~ NA_character_,

    # Estuaries / sea are excluded from the terrestrial comparison.
    code %in% c(522, 523) ~ NA_character_,

    # Grassland
    code == 231 ~ "Grassland",
    code %in% c(321, 322, 323, 324) ~ "Grassland",

    # Cropland
    code %in% c(
      211, 212, 213,
      221, 222, 223,
      241, 242, 243, 244
    ) ~ "Cropland",

    # Forest
    code %in% c(311, 312, 313) ~ "Forest",

    # Artificial / urban land
    code >= 111 & code <= 142 ~ "Urban",

    # Sparse/bare land, wetlands and inland water
    TRUE ~ "Other"
  )
}

# ============================================================
# 3. GBIF OCCURRENCE FILTERING
# ============================================================

occ_raw <- read_delim(
  GBIF_IN,
  delim = "\t",
  show_col_types = FALSE,
  guess_max = 100000
)

required_cols <- c(
  "decimalLatitude",
  "decimalLongitude",
  "occurrenceStatus",
  "basisOfRecord",
  "coordinateUncertaintyInMeters",
  "eventDate"
)

missing_cols <- setdiff(required_cols, names(occ_raw))

if (length(missing_cols) > 0) {
  stop(
    "Missing required GBIF columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

steps <- tibble(
  step = "Downloaded records",
  n = nrow(occ_raw)
)

add_count <- function(label, d) {
  steps <<- bind_rows(
    steps,
    tibble(step = label, n = nrow(d))
  )
}

occ <- occ_raw |>
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  )

add_count("Has coordinates", occ)

occ <- occ |>
  filter(occurrenceStatus == "PRESENT")

add_count("Occurrence status = PRESENT", occ)

occ <- occ |>
  filter(
    basisOfRecord %in%
      c("HUMAN_OBSERVATION", "OBSERVATION")
  )

add_count("Human observation", occ)

# A slightly stricter threshold than the GBIF download filter
# (which allowed uncertainty up to 1060 m).
occ <- occ |>
  filter(
    !is.na(coordinateUncertaintyInMeters),
    coordinateUncertaintyInMeters <= MAX_COORD_UNCERTAINTY
  )

add_count(
  paste0(
    "Coordinate uncertainty <= ",
    MAX_COORD_UNCERTAINTY,
    " m"
  ),
  occ
)

# Extract calendar year/month from the ISO eventDate field.
occ <- occ |>
  mutate(
    event_date = as.Date(substr(eventDate, 1, 10)),
    analysis_year = as.integer(format(event_date, "%Y")),
    analysis_month = as.integer(format(event_date, "%m"))
  ) |>
  filter(
    !is.na(analysis_year),
    analysis_year >= YEAR_MIN,
    analysis_year <= YEAR_MAX
  )

add_count(
  paste0("Years ", YEAR_MIN, "-", YEAR_MAX),
  occ
)

occ <- occ |>
  filter(analysis_month %in% FLIGHT_MONTHS)

add_count("Flight period (May-August)", occ)

occ <- occ |>
  distinct(
    decimalLatitude,
    decimalLongitude,
    event_date,
    .keep_all = TRUE
  )

add_count("Duplicate records removed", occ)

cat("\n==================== GBIF RECORD FILTERING ====================\n")
print(as.data.frame(steps), row.names = FALSE)

cat(
  "\nFinal year range:",
  min(occ$analysis_year),
  "-",
  max(occ$analysis_year),
  "\n"
)

write_tsv(
  steps,
  file.path(TAB, "gbif_record_filtering.tsv")
)

# ============================================================
# 4. PROJECT OCCURRENCES AND SPATIALLY THIN
# ============================================================

# Convert WGS84 GBIF coordinates to SWEREF99 TM so the same
# study-region bounds and 1-km grid can be used.

pts_ll <- vect(
  occ,
  geom = c("decimalLongitude", "decimalLatitude"),
  crs = "EPSG:4326"
)

pts_sweref <- project(pts_ll, SITE_CRS)
xy <- crds(pts_sweref)

occ$easting <- xy[, 1]
occ$northing <- xy[, 2]

occ <- occ |>
  mutate(
    in_study_region =
      easting >= STUDY_XMIN &
      easting <= STUDY_XMAX &
      northing >= STUDY_YMIN &
      northing <= STUDY_YMAX
  )

cat(
  "\nRecords in study region before thinning:",
  sum(occ$in_study_region),
  "of",
  nrow(occ),
  "\n"
)

# Thin to one record per 1-km SWEREF99 TM grid cell.
# This prevents repeatedly sampled localities from dominating
# the habitat proportions.

thin_cell <- paste(
  floor(occ$easting / THIN_GRID_M),
  floor(occ$northing / THIN_GRID_M)
)

keep <- ave(
  seq_along(thin_cell),
  thin_cell,
  FUN = function(k) k == sample(k, 1)
) == 1

occ_thin <- occ[keep, ]

cat(
  "After thinning to one record per",
  THIN_GRID_M / 1000,
  "km cell:",
  nrow(occ_thin),
  "records (",
  sum(occ_thin$in_study_region),
  "in study region)\n"
)

# ============================================================
# 5. CORINE LAND COVER
# ============================================================

if (!file.exists(CORINE_TIF)) {
  stop("CORINE raster not found: ", CORINE_TIF)
}

corine <- rast(CORINE_TIF)

clc_lookup <- read_tsv(
  CLC_LOOKUP,
  show_col_types = FALSE
) |>
  select(corine1, clc_code, clc_name) |>
  distinct()

index_to_group <- clc_lookup |>
  mutate(
    group = clc_to_group(clc_code)
  ) |>
  select(
    corine1,
    clc_code,
    clc_name,
    group
  )

# Crop CORINE to the same study-region extent used in the
# original analysis.
study_ext <- ext(
  STUDY_XMIN,
  STUDY_XMAX,
  STUDY_YMIN,
  STUDY_YMAX
)

study_poly <- as.polygons(
  study_ext,
  crs = SITE_CRS
) |>
  project(crs(corine))

corine_study <- crop(
  corine,
  study_poly,
  mask = TRUE
)

cat(
  "\nCropped CORINE:",
  paste(dim(corine_study)[1:2], collapse = " x "),
  "cells\n"
)

extract_group <- function(v, raster) {

  idx <- as.integer(
    terra::extract(raster, v)[, 2]
  )

  index_to_group$group[
    match(idx, index_to_group$corine1)
  ]
}

# ============================================================
# 6. GBIF PRESENCES
# ============================================================

pres_dat <- occ_thin |>
  filter(in_study_region)

pres_v <- vect(
  pres_dat,
  geom = c("easting", "northing"),
  crs = SITE_CRS
) |>
  project(crs(corine))

pres_group <- extract_group(
  pres_v,
  corine
)

# ============================================================
# 7. BACKGROUND A -- EXACT LANDSCAPE AVAILABILITY
# ============================================================

# Calculate habitat availability from every CORINE cell in the
# study region rather than from a random sample of cells.
# Because CORINE cells have equal area, cell proportions are
# equivalent to area proportions.

# CORINE is a categorical raster. terra::freq() therefore returns
# the active category labels rather than the underlying integer cell IDs.
# Recover the integer raster index from the raster's category table.

land_freq <- as.data.frame(
  terra::freq(corine_study)
)

if (!all(c("value", "count") %in% names(land_freq))) {
  stop("Unexpected output from terra::freq().")
}

corine_levels <- terra::levels(corine_study)[[1]]

if (is.null(corine_levels) || ncol(corine_levels) < 2) {
  stop("Could not retrieve CORINE raster category levels.")
}

id_col <- names(corine_levels)[1]
label_col <- names(corine_levels)[2]

land_freq$corine1 <- corine_levels[[id_col]][
  match(
    as.character(land_freq$value),
    as.character(corine_levels[[label_col]])
  )
]

if (any(is.na(land_freq$corine1))) {
  print(
    land_freq |>
      filter(is.na(corine1))
  )
  stop("Some CORINE category labels could not be mapped to raster indices.")
}

landscape_counts <- land_freq |>
  transmute(
    corine1 = as.integer(corine1),
    n_cells = count
  ) |>
  left_join(
    index_to_group,
    by = "corine1"
  ) |>
  filter(!is.na(group)) |>
  group_by(group) |>
  summarise(
    n_cells = sum(n_cells),
    .groups = "drop"
  )

landscape_availability <- tibble(
  habitat = GROUP_ORDER
) |>
  left_join(
    landscape_counts,
    by = c("habitat" = "group")
  ) |>
  mutate(
    n_cells = coalesce(n_cells, 0),
    proportion = n_cells / sum(n_cells),
    percent = 100 * proportion
  )

cat(
  "\n==================== LANDSCAPE AVAILABILITY ====================\n"
)

print(
  as.data.frame(
    landscape_availability |>
      mutate(percent = round(percent, 2))
  ),
  row.names = FALSE
)

# ============================================================
# 8. BACKGROUND B -- SEBMS MONITORING SITES
# ============================================================

# Monitoring sites provide a second background that partially
# accounts for accessibility / observer-location bias.

sites <- read_tsv(
  SITES_IN,
  show_col_types = FALSE
)

required_site_cols <- c(
  "in_study_region",
  "clc_code_extracted"
)

missing_site_cols <- setdiff(
  required_site_cols,
  names(sites)
)

if (length(missing_site_cols) > 0) {
  stop(
    "Missing required site columns: ",
    paste(missing_site_cols, collapse = ", ")
  )
}

bg_sites_group <- sites |>
  filter(in_study_region) |>
  mutate(
    group = clc_to_group(clc_code_extracted)
  ) |>
  pull(group)

cat(
  "\n==================== ANALYSIS SAMPLE SIZES ====================\n"
)

cat(
  "GBIF occurrence records (study region, thinned, terrestrial):",
  sum(!is.na(pres_group)),
  "\n"
)

cat(
  "Terrestrial CORINE cells in study region:",
  format(
    sum(landscape_availability$n_cells),
    big.mark = ","
  ),
  "\n"
)

cat(
  "SeBMS monitoring sites (terrestrial):",
  sum(!is.na(bg_sites_group)),
  "\n"
)

# ============================================================
# 9. SELECTION RATIOS
# ============================================================

# Selection ratio =
# proportion of GBIF records in habitat /
# proportion of background locations in habitat.
#
# Values >1 indicate overrepresentation relative to that
# background; values <1 indicate underrepresentation.
#
# Opportunistic GBIF records can be affected by spatially
# heterogeneous recorder effort. Therefore interpretation
# focuses primarily on habitat classes showing consistent
# patterns against both backgrounds.

habitat_props <- function(g) {

  g <- g[!is.na(g)]

  tab <- table(
    factor(g, levels = GROUP_ORDER)
  )

  as.numeric(tab) / sum(tab)
}

# Landscape availability is known from the complete CORINE
# raster, so only occurrence records are bootstrapped.
selection_ratios_landscape <- function(
    pres,
    background_props,
    label
) {

  pres <- pres[!is.na(pres)]

  p_obs <- habitat_props(pres)

  ratio_obs <- p_obs / background_props

  boot <- replicate(
    N_BOOT,
    {
      p <- habitat_props(
        sample(
          pres,
          length(pres),
          replace = TRUE
        )
      )

      p / background_props
    }
  )

  tibble(
    background = label,
    habitat = GROUP_ORDER,
    n_records = as.numeric(
      table(
        factor(
          pres,
          levels = GROUP_ORDER
        )
      )
    ),
    pct_records = 100 * p_obs,
    pct_background = 100 * background_props,
    ratio = ratio_obs,
    ratio_low = apply(
      boot,
      1,
      quantile,
      0.025,
      na.rm = TRUE
    ),
    ratio_high = apply(
      boot,
      1,
      quantile,
      0.975,
      na.rm = TRUE
    )
  )
}

# For the monitoring-site comparison both the occurrence and
# monitoring-site samples are bootstrapped.
selection_ratios_monitoring <- function(
    pres,
    bg,
    label
) {

  pres <- pres[!is.na(pres)]
  bg <- bg[!is.na(bg)]

  p_obs <- habitat_props(pres)
  b_obs <- habitat_props(bg)

  ratio_obs <- p_obs / b_obs

  boot <- replicate(
    N_BOOT,
    {
      p <- habitat_props(
        sample(
          pres,
          length(pres),
          replace = TRUE
        )
      )

      b <- habitat_props(
        sample(
          bg,
          length(bg),
          replace = TRUE
        )
      )

      p / b
    }
  )

  tibble(
    background = label,
    habitat = GROUP_ORDER,
    n_records = as.numeric(
      table(
        factor(
          pres,
          levels = GROUP_ORDER
        )
      )
    ),
    pct_records = 100 * p_obs,
    pct_background = 100 * b_obs,
    ratio = ratio_obs,
    ratio_low = apply(
      boot,
      1,
      quantile,
      0.025,
      na.rm = TRUE
    ),
    ratio_high = apply(
      boot,
      1,
      quantile,
      0.975,
      na.rm = TRUE
    )
  )
}

ratios <- bind_rows(
  selection_ratios_landscape(
    pres_group,
    landscape_availability$proportion,
    "Landscape availability"
  ),
  selection_ratios_monitoring(
    pres_group,
    bg_sites_group,
    "Monitoring sites"
  )
)

write_tsv(
  ratios,
  file.path(
    TAB,
    "TABLE_S4_gbif_habitat_selection_numeric.tsv"
  )
)

cat(
  "\n==================== SELECTION RATIOS ====================\n"
)

print(
  as.data.frame(
    ratios |>
      mutate(
        across(
          where(is.numeric),
          \(x) round(x, 3)
        )
      )
  ),
  row.names = FALSE
)

# ============================================================
# 10. CONSISTENCY ACROSS BACKGROUNDS
# ============================================================

direction_from_ci <- function(lo, hi) {
  case_when(
    lo > 1 ~ "above 1",
    hi < 1 ~ "below 1",
    TRUE ~ "not distinguishable from 1"
  )
}

landscape_direction <- ratios |>
  filter(background == "Landscape availability") |>
  transmute(
    habitat,
    landscape_direction =
      direction_from_ci(ratio_low, ratio_high)
  )

monitor_direction <- ratios |>
  filter(background == "Monitoring sites") |>
  transmute(
    habitat,
    monitoring_direction =
      direction_from_ci(ratio_low, ratio_high)
  )

agreement <- left_join(
  landscape_direction,
  monitor_direction,
  by = "habitat"
) |>
  mutate(
    consistent =
      landscape_direction == monitoring_direction
  )

cat(
  "\n==================== DIRECTION ACROSS BACKGROUNDS ====================\n"
)

print(
  as.data.frame(agreement),
  row.names = FALSE
)

# ============================================================
# 11. TABLE S4
# ============================================================

fmt_ratio <- function(r, lo, hi) {
  sprintf(
    "%.2f (%.2f-%.2f)",
    r,
    lo,
    hi
  )
}

random_tab <- ratios |>
  filter(background == "Landscape availability") |>
  transmute(
    `Habitat class` = habitat,

    `Occurrence records, n (%)` =
      sprintf(
        "%d (%.1f)",
        n_records,
        pct_records
      ),

    `Landscape availability (%)` =
      sprintf("%.1f", pct_background),

    `Landscape selection ratio (95% CI)` =
      fmt_ratio(
        ratio,
        ratio_low,
        ratio_high
      )
  )

monitor_tab <- ratios |>
  filter(background == "Monitoring sites") |>
  transmute(
    `Habitat class` = habitat,

    `Monitoring sites (%)` =
      sprintf("%.1f", pct_background),

    `Monitoring-site selection ratio (95% CI)` =
      fmt_ratio(
        ratio,
        ratio_low,
        ratio_high
      )
  )

table_s4 <- left_join(
  random_tab,
  monitor_tab,
  by = "Habitat class"
)

write_tsv(
  table_s4,
  file.path(
    TAB,
    "TABLE_S4_gbif_habitat_selection.tsv"
  )
)

cat(
  "\n==================== TABLE S4 ====================\n"
)

print(
  as.data.frame(table_s4),
  right = FALSE,
  row.names = FALSE
)

# ============================================================
# 12. OUTPUT SUMMARY
# ============================================================

cat(
  "\n==================== OUTPUTS ====================\n"
)

cat(
  "Table S4:\n",
  "  TABLE_S4_gbif_habitat_selection.tsv\n",
  "  TABLE_S4_gbif_habitat_selection_numeric.tsv\n",
  sep = ""
)

cat(
  "\nFiltering record:\n",
  "  gbif_record_filtering.tsv\n",
  sep = ""
)

cat(
  "\nTables -> ",
  TAB,
  "\n",
  sep = ""
)
