# ============================================================
# Generate annual HILDA+ habitat suitability maps
#
# Inputs:
#   01-data/HILDA/preprocessed_wgs84/
#   01-data/HILDA/hilda_final_clip_20km_south.tif
#   01-data/boundaries/sweden_shape.shp
#
# Outputs:
#   03-model_input/suit_maps/hilda_suitability_<year>.tif
#   03-model_input/suit_maps/hilda_suitability_<year>.txt
#   03-model_input/suit_maps/hilda_suitability_<year>_dims.txt
#
#   02-landscapes/derived/class_areas_by_year.tsv
#   02-landscapes/derived/suitability_value_counts.tsv
#
# Usage:
#   Rscript 00-scripts/03f_make_suitability_maps.R
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
})

# ============================================================
# 1. PATHS
# ============================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Run this script using Rscript.")
}

SCRIPT_DIR <- dirname(
  normalizePath(sub("^--file=", "", file_arg[1]))
)

REPO_ROOT <- normalizePath(
  file.path(SCRIPT_DIR, "..")
)

HILDA_DIR <- file.path(
  REPO_ROOT,
  "01-data",
  "HILDA",
  "preprocessed_wgs84"
)

TEMPLATE <- file.path(
  REPO_ROOT,
  "01-data",
  "HILDA",
  "hilda_final_clip_20km_south.tif"
)

SWEDEN_SHP <- file.path(
  REPO_ROOT,
  "01-data",
  "boundaries",
  "sweden_shape.shp"
)

OUT_DIR <- file.path(
  REPO_ROOT,
  "03-model_input",
  "suit_maps"
)

DER_DIR <- file.path(
  REPO_ROOT,
  "02-landscapes",
  "derived"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DER_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. HABITAT SUITABILITY
# ============================================================

RECLASS <- tibble::tribble(
  ~code, ~hilda_label,                   ~model_class,  ~suitability,
   0L,   "Ocean / no data",              "Unsuitable",  0.0,
  11L,   "Urban",                        "Urban",        0.1,
  22L,   "Cropland",                     "Cropland",     0.5,
  33L,   "Pasture / rangeland",          "Grassland",    1.0,
  44L,   "Forest",                       "Forest",       0.3,
  55L,   "Unmanaged grass / shrubland",  "Grassland",    1.0,
  66L,   "Sparse / no vegetation",       "Unsuitable",  0.0,
  77L,   "Water",                        "Unsuitable",   0.0
)

GRASSLAND_VALUE <- 1.0

# ============================================================
# 3. INPUTS
# ============================================================

if (!dir.exists(HILDA_DIR)) {
  stop("HILDA directory not found: ", HILDA_DIR)
}

if (!file.exists(TEMPLATE)) {
  stop("Template raster not found: ", TEMPLATE)
}

if (!file.exists(SWEDEN_SHP)) {
  stop("Sweden boundary not found: ", SWEDEN_SHP)
}

template <- rast(TEMPLATE)

sweden <- vect(SWEDEN_SHP)

if (!same.crs(template, sweden)) {
  sweden <- project(
    sweden,
    crs(template)
  )
}

EXPECTED_NCELL <- 94180L

if (ncell(template) != EXPECTED_NCELL) {
  stop(
    "Template has ",
    ncell(template),
    " cells; expected ",
    EXPECTED_NCELL
  )
}

cell_km2 <- prod(res(template)) / 1e6

cat("\n==================== SIMULATION GRID ====================\n")

cat(
  "Dimensions :",
  nrow(template),
  "x",
  ncol(template),
  "=",
  ncell(template),
  "cells\n"
)

cat(
  "Resolution :",
  round(res(template)[1], 6),
  "x",
  round(res(template)[2], 6),
  "m\n"
)

cat(
  "Cell area  :",
  round(cell_km2, 6),
  "km2\n"
)

cat(
  "Extent     :",
  xmin(template),
  xmax(template),
  ymin(template),
  ymax(template),
  "\n"
)

# ============================================================
# 4. FIND ANNUAL HILDA+ FILES
# ============================================================

HILDA_PATTERN <-
  "^hilda_plus_([0-9]{4})_states_GLOB-v1-0_wgs84-nn\\.tif$"

hilda_files <- list.files(
  HILDA_DIR,
  pattern = HILDA_PATTERN,
  full.names = TRUE
)

years <- as.integer(
  sub(
    HILDA_PATTERN,
    "\\1",
    basename(hilda_files)
  )
)

ord <- order(years)

years <- years[ord]
hilda_files <- hilda_files[ord]

EXPECTED_YEARS <- 1899:2019

missing_years <- setdiff(
  EXPECTED_YEARS,
  years
)

if (length(missing_years) > 0) {
  stop(
    "Missing HILDA+ years: ",
    paste(missing_years, collapse = ", ")
  )
}

keep <- years %in% EXPECTED_YEARS

years <- years[keep]
hilda_files <- hilda_files[keep]

cat(
  "\nHILDA+ years:",
  min(years),
  "-",
  max(years),
  "(",
  length(years),
  "rasters)\n"
)

# ============================================================
# 5. BUILD ANNUAL MAPS
# ============================================================

area_log <- list()
value_log <- list()

cat("\n==================== BUILDING MAPS ====================\n")

for (i in seq_along(years)) {

  yr <- years[i]

  r_raw <- rast(
    hilda_files[i]
  )

  # Project the preprocessed WGS84 HILDA raster onto the
  # simulation template using nearest-neighbour interpolation.
  r_match <- project(
    r_raw,
    template,
    method = "near"
  )

  # Mask to Sweden using the same sequence used to construct
  # the original simulation landscapes.
  r_crop <- crop(
    r_match,
    sweden
  )

  r_swe <- mask(
    r_crop,
    sweden
  )

  r_swe <- extend(
    r_swe,
    ext(template)
  )

  r_swe <- resample(
    r_swe,
    template,
    method = "near"
  )

  r_swe[is.na(r_swe)] <- 0

  vr <- values(
    r_swe,
    mat = FALSE
  )

  # Validate land-cover classes.
  unexpected_codes <- setdiff(
    unique(as.integer(vr)),
    RECLASS$code
  )

  if (length(unexpected_codes) > 0) {
    stop(
      "Unexpected HILDA+ classes in ",
      yr,
      ": ",
      paste(unexpected_codes, collapse = ", ")
    )
  }

  # Original land-cover composition.
  tab_raw <- table(vr)

  area_log[[as.character(yr)]] <- tibble(
    year = yr,
    code = as.integer(names(tab_raw)),
    cells = as.integer(tab_raw)
  )

  # Reclassify to relative suitability.
  suit <- subst(
    r_swe,
    from = RECLASS$code,
    to = RECLASS$suitability
  )

  suit[is.na(suit)] <- 0

  vs <- values(
    suit,
    mat = FALSE
  )

  if (length(vs) != EXPECTED_NCELL) {
    stop(
      "Incorrect number of cells in ",
      yr
    )
  }

  unexpected_suit <- setdiff(
    unique(vs),
    unique(RECLASS$suitability)
  )

  if (length(unexpected_suit) > 0) {
    stop(
      "Unexpected suitability values in ",
      yr,
      ": ",
      paste(unexpected_suit, collapse = ", ")
    )
  }

  tab_suit <- table(vs)

  value_log[[as.character(yr)]] <- tibble(
    year = yr,
    suitability = as.numeric(names(tab_suit)),
    cells = as.integer(tab_suit)
  )

  # ----------------------------------------------------------
  # Write simulation inputs
  # ----------------------------------------------------------

  writeRaster(
    suit,
    file.path(
      OUT_DIR,
      sprintf(
        "hilda_suitability_%d.tif",
        yr
      )
    ),
    overwrite = TRUE
  )

  # terra returns raster values in row-major order, matching
  # matrix(..., byrow = TRUE) when read by SLiM.
  writeLines(
    format(
      vs,
      scientific = FALSE,
      trim = TRUE
    ),
    file.path(
      OUT_DIR,
      sprintf(
        "hilda_suitability_%d.txt",
        yr
      )
    )
  )

  writeLines(
    c(
      paste0("year\t", yr),
      paste0("nrow\t", nrow(suit)),
      paste0("ncol\t", ncol(suit)),
      paste0("xmin\t", xmin(suit)),
      paste0("xmax\t", xmax(suit)),
      paste0("ymin\t", ymin(suit)),
      paste0("ymax\t", ymax(suit)),
      paste0("xres\t", res(suit)[1]),
      paste0("yres\t", res(suit)[2])
    ),
    file.path(
      OUT_DIR,
      sprintf(
        "hilda_suitability_%d_dims.txt",
        yr
      )
    )
  )

  grass_cells <- sum(
    vs == GRASSLAND_VALUE
  )

  cat(
    sprintf(
      "%d  grassland: %5d cells  (%7.1f km2)\n",
      yr,
      grass_cells,
      grass_cells * cell_km2
    )
  )
}

# ============================================================
# 6. LAND-COVER SUMMARY
# ============================================================

areas <- bind_rows(
  area_log
) |>
  left_join(
    RECLASS |>
      select(
        code,
        hilda_label,
        model_class,
        suitability
      ),
    by = "code"
  ) |>
  mutate(
    area_km2 = cells * cell_km2
  ) |>
  arrange(
    year,
    code
  )

write_tsv(
  areas,
  file.path(
    DER_DIR,
    "class_areas_by_year.tsv"
  )
)

# ============================================================
# 7. SUITABILITY SUMMARY
# ============================================================

vals <- bind_rows(
  value_log
) |>
  mutate(
    area_km2 = cells * cell_km2
  ) |>
  arrange(
    year,
    suitability
  )

write_tsv(
  vals,
  file.path(
    DER_DIR,
    "suitability_value_counts.tsv"
  )
)

grass <- vals |>
  filter(
    suitability == GRASSLAND_VALUE
  )

g_first <- grass |>
  filter(
    year == min(year)
  )

g_last <- grass |>
  filter(
    year == max(year)
  )

cat(
  "\n==================== GRASSLAND TRAJECTORY ====================\n"
)

cat(
  sprintf(
    "%d : %.1f km2\n",
    g_first$year,
    g_first$area_km2
  )
)

cat(
  sprintf(
    "%d : %.1f km2\n",
    g_last$year,
    g_last$area_km2
  )
)

cat(
  sprintf(
    "Change: %+.1f%%\n",
    100 *
      (g_last$area_km2 - g_first$area_km2) /
      g_first$area_km2
  )
)

cat(
  "\nMaps -> ",
  OUT_DIR,
  "\n",
  sep = ""
)

cat(
  "Tables -> ",
  DER_DIR,
  "\n",
  sep = ""
)
