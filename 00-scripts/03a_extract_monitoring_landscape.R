
# ============================================================
# Extract habitat and landscape composition at monitoring sites
#
# Inputs:
#   01-data/monitoring/presence_absence_master.tsv
#   01-data/monitoring/corine_index_to_clc.tsv
#   01-data/CORINE/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif
#   01-data/HILDA/hildap_vGLOB-1.0_geotiff_eckert4/states/hilda_plus_2019_states_GLOB-v20-11-29.tif
#
# Outputs:
#   02-landscapes/derived/sites_with_classes.tsv
#   02-landscapes/derived/site_buffer_composition_long.tsv
#
# Usage:
#   Rscript 00-scripts/03a_extract_monitoring_landscape.R
# ============================================================

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(tidyr); library(readr)
})

# ==========================================================================
# 1. PATHS
# ==========================================================================
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Run this script using Rscript.")
}

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))

MONITORING_DIR <- file.path(REPO_ROOT, "01-data", "monitoring")
OUT            <- file.path(REPO_ROOT, "02-landscapes", "derived")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SITES_IN   <- file.path(MONITORING_DIR, "presence_absence_master.tsv")
CLC_LOOKUP <- file.path(MONITORING_DIR, "corine_index_to_clc.tsv")

CORINE_TIF <- file.path(
  REPO_ROOT,
  "01-data", "CORINE",
  "u2018_clc2018_v2020_20u1_raster100m",
  "DATA",
  "U2018_CLC2018_V2020_20u1.tif"
)

HILDA_TIF <- file.path(
  REPO_ROOT,
  "01-data", "HILDA",
  "hildap_vGLOB-1.0_geotiff_eckert4",
  "states",
  "hilda_plus_2019_states_GLOB-v20-11-29.tif"
)

SITE_CRS   <- "EPSG:3006"        # SWEREF99 TM; lon = easting, lat = northing
BUF_CORINE <- c(500, 1000, 2000, 5000, 10000)
BUF_HILDA  <- c(2000, 5000, 10000, 20000)

# Sites processed per extraction chunk. Lower this if memory is tight; it
# changes only speed and peak memory, never the result.
CHUNK_SITES <- 150

# Study region: the simulated landscape extent, in EPSG:3006.
# Used only to flag sites, not to drop them.
STUDY_YMIN <- 6112797; STUDY_YMAX <- 6418978
STUDY_XMIN <-  295217; STUDY_XMAX <-  672595

HILDA_LABELS <- c(`0` = "ocean / no data", `11` = "urban", `22` = "cropland",
                  `33` = "pasture / rangeland", `44` = "forest",
                  `55` = "unmanaged grass / shrubland",
                  `66` = "sparse / no vegetation", `77` = "water")

# ==========================================================================
# 2. SITES
# ==========================================================================
sites <- read_tsv(SITES_IN, show_col_types = FALSE)
stopifnot(all(c("sit_uid", "lat", "lon", "presence") %in% names(sites)))

# --- sit_uid must be unique. Everything downstream joins on it, so a repeated
# --- value silently duplicates rows in the buffer models and turns one site
# --- into a pseudo-replicate. Inspect before continuing rather than guessing
# --- which copy to keep.
dup_uid <- sites |> count(sit_uid) |> filter(n > 1)

if (nrow(dup_uid) > 0) {
  cat("!!! duplicated sit_uid values:", nrow(dup_uid), "\n")
  print(as.data.frame(
    sites |> filter(sit_uid %in% dup_uid$sit_uid) |> arrange(sit_uid)
  ))
  stop("sit_uid is not unique. Resolve these rows before extracting.")
}

cat("sit_uid unique:", nrow(sites), "sites\n")

# lon = easting, lat = northing, despite the column names
pts <- vect(as.data.frame(sites), geom = c("lon", "lat"), crs = SITE_CRS)

sites <- sites |>
  mutate(in_study_region = lon >= STUDY_XMIN & lon <= STUDY_XMAX &
           lat >= STUDY_YMIN & lat <= STUDY_YMAX)

cat("sites:", nrow(sites),
    " presence:", sum(sites$presence),
    " in study region:", sum(sites$in_study_region), "\n\n")

# ==========================================================================
# 3. RASTERS
# ==========================================================================
stopifnot(file.exists(CORINE_TIF), file.exists(HILDA_TIF))
corine <- rast(CORINE_TIF)
hilda  <- rast(HILDA_TIF)

cat("CORINE CRS:", crs(corine, describe = TRUE)$name,
    " res:", paste(res(corine), collapse = " x "), "\n")
cat("HILDA  CRS:", crs(hilda,  describe = TRUE)$name,
    " res:", paste(res(hilda),  collapse = " x "), "\n\n")

# The CORINE raster carries a category table, so extract() returns a factor
# whose integer codes are the raster index values 1-44 (column `corine1`).
clc_lookup <- read_tsv(CLC_LOOKUP, show_col_types = FALSE) |>
  select(corine1, clc_code, clc_name) |>
  distinct()

# ==========================================================================
# 4. POINT-LEVEL CLASSES
# ==========================================================================
pt_clc   <- terra::extract(corine, project(pts, crs(corine)))[, 2]
pt_hilda <- terra::extract(hilda,  project(pts, crs(hilda)))[, 2]

sites <- sites |>
  mutate(clc_index_extracted = as.integer(pt_clc),
         hilda_code          = as.integer(pt_hilda),
         hilda_class         = unname(HILDA_LABELS[as.character(hilda_code)])) |>
  left_join(clc_lookup |>
              select(corine1,
                     clc_code_extracted = clc_code,
                     clc_name_extracted = clc_name),
            by = c("clc_index_extracted" = "corine1"))

# --- check 1: does my extraction reproduce the CORINE class already there?
if ("clc_code" %in% names(sites)) {
  agree <- sum(as.integer(sites$clc_code_extracted) == as.integer(sites$clc_code),
               na.rm = TRUE)
  cat("--- CORINE cross-check ---\n")
  cat("  agreement with existing clc_code:", agree, "/", nrow(sites),
      sprintf("(%.1f%%)\n", 100 * agree / nrow(sites)))
  cat("\n")
}

occ_table <- function(d) {
  d |>
    count(hilda_code, hilda_class, presence) |>
    pivot_wider(names_from = presence, values_from = n,
                names_prefix = "presence_", values_fill = 0) |>
    mutate(n_sites   = presence_0 + presence_1,
           occupancy = round(presence_1 / n_sites, 3)) |>
    arrange(desc(n_sites))
}

cat("--- HILDA class at site, by presence (all sites) ---\n")
print(as.data.frame(occ_table(sites)))
cat("\n--- HILDA class at site, by presence (study region only) ---\n")
print(as.data.frame(occ_table(sites |> filter(in_study_region))))
cat("\n")

# ==========================================================================
# 5. BUFFER COMPOSITION
#    One row per site per class, with the proportion of the buffer.
#    as.integer() on the CORINE factor recovers the 1-44 raster index, which
#    is what the point-level extraction used and what `corine1` holds.
#
#    Sites are processed in chunks and each chunk is reduced to class counts
#    before the next is read, so peak memory does not scale with the number
#    of sites. Results are identical to extracting everything at once.
# ==========================================================================
buffer_composition <- function(points, raster, radius_m, source_label,
                               chunk_size = CHUNK_SITES) {
  n_pts  <- nrow(points)
  starts <- seq(1, n_pts, by = chunk_size)
  t0     <- Sys.time()

  chunks <- lapply(seq_along(starts), function(k) {
    i0 <- starts[k]
    i1 <- min(i0 + chunk_size - 1, n_pts)

    sub  <- points[i0:i1, ]
    buf  <- buffer(sub, width = radius_m)
    bufp <- project(buf, crs(raster))

    ex <- terra::extract(raster, bufp)
    names(ex) <- c("ID", "value")
    ex$value <- as.integer(ex$value)   # factor (CORINE) or integer (HILDA)

    out <- ex |>
      filter(!is.na(value)) |>
      count(ID, value, name = "cells") |>
      mutate(sit_uid = sub$sit_uid[ID]) |>
      select(sit_uid, value, cells)

    rm(ex, buf, bufp); gc(verbose = FALSE)

    cat(sprintf("\r  %s %g m: sites %d-%d of %d (%.0f s)   ",
                source_label, radius_m, i0, i1, n_pts,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush.console()

    out
  })

  cat("\n")

  bind_rows(chunks) |>
    group_by(sit_uid, value) |>
    summarise(cells = sum(cells), .groups = "drop") |>
    group_by(sit_uid) |>
    mutate(prop        = cells / sum(cells),
           total_cells = sum(cells)) |>
    ungroup() |>
    mutate(radius_m = radius_m,
           source   = source_label) |>
    select(sit_uid, source, radius_m, value, cells, total_cells, prop)
}

comp <- bind_rows(
  lapply(BUF_CORINE, function(r) {
    cat("CORINE buffer", r, "m ...\n")
    buffer_composition(pts, corine, r, "corine")
  }),
  lapply(BUF_HILDA, function(r) {
    cat("HILDA  buffer", r, "m ...\n")
    buffer_composition(pts, hilda, r, "hilda")
  })
)

# translate raw raster values into named classes
comp <- comp |>
  left_join(clc_lookup, by = c("value" = "corine1")) |>
  mutate(clc_code    = if_else(source == "corine", clc_code, NA_real_),
         clc_name    = if_else(source == "corine", clc_name, NA_character_),
         hilda_class = if_else(source == "hilda",
                               unname(HILDA_LABELS[as.character(value)]),
                               NA_character_))

write_tsv(comp, file.path(OUT, "site_buffer_composition_long.tsv"))

# --- sanity check 1: proportions sum to 1 per site per radius
chk <- comp |>
  group_by(sit_uid, source, radius_m) |>
  summarise(total_prop = sum(prop), .groups = "drop")
cat("\nproportion sums: min", round(min(chk$total_prop), 6),
    " max", round(max(chk$total_prop), 6), "\n")

# --- sanity check 2: how many sites came back per source x radius, and the
# --- median buffer size in cells. Sites missing at large radii are the ones
# --- whose buffer is almost entirely sea; the analysis script handles those,
# --- but the counts should be inspected here first.
cat("\n--- sites and buffer size per source x radius ---\n")
print(as.data.frame(
  comp |>
    group_by(source, radius_m) |>
    summarise(n_sites = n_distinct(sit_uid),
              median_cells = median(total_cells),
              .groups = "drop")
))

## Check where corine sites go to in Hilda classification
if ("habitat_hilda" %in% names(sites)) {
  cat("\n--- CORINE-derived group vs HILDA class at site ---\n")
  print(
    sites |>
      count(habitat_hilda, hilda_class) |>
      pivot_wider(names_from = hilda_class, values_from = n, values_fill = 0) |>
      as.data.frame()
  )
}

# ==========================================================================
# 6. WRITE SITE TABLE
# ==========================================================================
write_tsv(sites, file.path(OUT, "sites_with_classes.tsv"))

cat("\nwrote:\n")
cat("  ", file.path(OUT, "sites_with_classes.tsv"), "\n", sep = "")
cat("  ", file.path(OUT, "site_buffer_composition_long.tsv"),
    " (", nrow(comp), " rows)\n", sep = "")
