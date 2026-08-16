# ============================================================
# Summarize historical habitat change, 1899-2019
#
# Input:
#   02-landscapes/derived/class_areas_by_year.tsv
#
# Outputs:
#   02-landscapes/figures/FIG_5_habitat_change_1899_2019.*
#   02-landscapes/tables/TABLE_S5_land_cover_change.tsv
#   02-landscapes/tables/TABLE_S5_land_cover_change_numeric.tsv
#   02-landscapes/tables/habitat_change_summary.tsv
#
# Usage:
#   Rscript 00-scripts/03g_habitat_change_summary.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
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

AREAS_IN <- file.path(
  REPO_ROOT,
  "02-landscapes",
  "derived",
  "class_areas_by_year.tsv"
)

FIG_DIR <- file.path(
  REPO_ROOT,
  "02-landscapes",
  "figures"
)

TAB_DIR <- file.path(
  REPO_ROOT,
  "02-landscapes",
  "tables"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TAB_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(AREAS_IN)) {
  stop(
    "Input not found: ",
    AREAS_IN,
    "\nRun 03f_make_suitability_maps.R first."
  )
}

# ============================================================
# 2. READ LAND-COVER AREAS
# ============================================================

areas <- read_tsv(
  AREAS_IN,
  show_col_types = FALSE
)

required_cols <- c(
  "year",
  "code",
  "cells",
  "hilda_label",
  "model_class",
  "suitability",
  "area_km2"
)

missing_cols <- setdiff(
  required_cols,
  names(areas)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing columns in class_areas_by_year.tsv: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (!all(1899:2019 %in% areas$year)) {
  stop("Expected annual data from 1899 through 2019.")
}

# Confirm one consistent cell area.
cell_area <- areas |>
  filter(cells > 0) |>
  transmute(
    cell_km2 = area_km2 / cells
  ) |>
  pull(cell_km2)

cell_km2 <- median(cell_area)

if (max(abs(cell_area - cell_km2)) > 1e-6) {
  stop("Inconsistent raster-cell areas detected.")
}

cat(
  "\nCell area:",
  round(cell_km2, 6),
  "km2\n"
)

# ============================================================
# 3. WEIGHTED SUITABLE HABITAT
# ============================================================

# Each cell contributes its area multiplied by its assigned
# relative suitability value.

weighted <- areas |>
  group_by(year) |>
  summarise(
    weighted_km2 =
      sum(area_km2 * suitability),
    .groups = "drop"
  )

# ============================================================
# 4. GRASSLAND AREA
# ============================================================

# HILDA+ classes 33 and 55 are both assigned suitability 1.0
# and treated as the model's highest-suitability grassland.

grass <- areas |>
  filter(
    code %in% c(33L, 55L)
  ) |>
  mutate(
    class = factor(
      code,
      levels = c(55L, 33L),
      labels = c(
        "Unmanaged grass / shrubland",
        "Pasture / rangeland"
      )
    )
  )

grass_total <- grass |>
  group_by(year) |>
  summarise(
    grassland_km2 = sum(area_km2),
    .groups = "drop"
  )

# ============================================================
# 5. SUMMARY OF HISTORICAL CHANGE
# ============================================================

weighted_1899 <- weighted |>
  filter(year == 1899) |>
  pull(weighted_km2)

weighted_2019 <- weighted |>
  filter(year == 2019) |>
  pull(weighted_km2)

grass_1899 <- grass_total |>
  filter(year == 1899) |>
  pull(grassland_km2)

grass_2019 <- grass_total |>
  filter(year == 2019) |>
  pull(grassland_km2)

summary_change <- tibble(
  metric = c(
    "Weighted suitable habitat",
    "Total grassland"
  ),
  area_1899_km2 = c(
    weighted_1899,
    grass_1899
  ),
  area_2019_km2 = c(
    weighted_2019,
    grass_2019
  )
) |>
  mutate(
    change_km2 =
      area_2019_km2 - area_1899_km2,

    change_percent =
      100 *
      (
        area_2019_km2 /
          area_1899_km2 -
          1
      )
  )

write_tsv(
  summary_change,
  file.path(
    TAB_DIR,
    "habitat_change_summary.tsv"
  )
)

# ============================================================
# 6. TABLE S5 -- LAND-COVER CHANGE
# ============================================================

# Percentages are calculated relative to terrestrial cells
# with suitability > 0, matching the thesis table.

land_totals <- areas |>
  filter(
    suitability > 0,
    year %in% c(1899, 2019)
  ) |>
  group_by(year) |>
  summarise(
    land_km2 = sum(area_km2),
    .groups = "drop"
  )

table_numeric <- areas |>
  filter(
    year %in% c(1899, 2019),
    code %in% c(55L, 33L, 22L, 44L, 11L)
  ) |>
  select(
    year,
    code,
    hilda_label,
    model_class,
    suitability,
    area_km2
  ) |>
  left_join(
    land_totals,
    by = "year"
  ) |>
  mutate(
    pct_land =
      100 * area_km2 / land_km2
  ) |>
  select(
    -land_km2
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(
      area_km2,
      pct_land
    )
  ) |>
  mutate(
    change_percent =
      100 *
      (
        area_km2_2019 /
          area_km2_1899 -
          1
      )
  ) |>
  arrange(
    match(
      code,
      c(55L, 33L, 22L, 44L, 11L)
    )
  )

write_tsv(
  table_numeric,
  file.path(
    TAB_DIR,
    "TABLE_S5_land_cover_change_numeric.tsv"
  )
)

fmt_change <- function(x) {
  sprintf("%+.1f%%", x)
}

table_formatted <- table_numeric |>
  transmute(
    `HILDA+ class (code)` =
      paste0(
        hilda_label,
        " (",
        code,
        ")"
      ),

    `Model class` =
      model_class,

    `Suitability value` =
      sprintf(
        "%.1f",
        suitability
      ),

    `1899 km2` =
      format(
        round(area_km2_1899),
        big.mark = ",",
        scientific = FALSE
      ),

    `% land 1899` =
      sprintf(
        "%.2f",
        pct_land_1899
      ),

    `2019 km2` =
      format(
        round(area_km2_2019),
        big.mark = ",",
        scientific = FALSE
      ),

    `% land 2019` =
      sprintf(
        "%.2f",
        pct_land_2019
      ),

    `Change` =
      fmt_change(
        change_percent
      )
  )

table_formatted <- bind_rows(
  table_formatted,
  tibble(
    `HILDA+ class (code)` =
      "Water, bare ground, outside area (0, 66, 77)",

    `Model class` =
      "Unsuitable",

    `Suitability value` =
      "0.0",

    `1899 km2` = "-",
    `% land 1899` = "-",
    `2019 km2` = "-",
    `% land 2019` = "-",
    `Change` = "-"
  )
)

write_tsv(
  table_formatted,
  file.path(
    TAB_DIR,
    "TABLE_S5_land_cover_change.tsv"
  )
)

# ============================================================
# 7. FIGURE 5
# ============================================================

# Colours retained from the final landscape figure.
COL_WEIGHTED  <- "#1B5E20"
COL_PASTURE   <- "#3E8E6B"
COL_UNMANAGED <- "#9CCC8F"

x_breaks <- seq(
  1900,
  2020,
  by = 20
)

# Panel A: weighted suitable habitat.
# The y-axis is intentionally fitted to the data rather than
# starting at zero, as stated in the figure caption.

pA <- ggplot(
  weighted,
  aes(
    x = year,
    y = weighted_km2
  )
) +
  geom_line(
    linewidth = 0.9,
    colour = COL_WEIGHTED
  ) +
  scale_x_continuous(
    breaks = x_breaks
  ) +
  labs(
    x = "Year",
    y = expression(
      "Available suitable habitat area (weighted km"^2*")"
    )
  ) +
  theme_classic(
    base_size = 12
  )

# Panel B: grassland classes shown as a stack so that the upper
# boundary represents total grassland area.

pB <- ggplot(
  grass,
  aes(
    x = year,
    y = area_km2,
    fill = class
  )
) +
  geom_area() +
  scale_fill_manual(
    values = c(
      "Unmanaged grass / shrubland" =
        COL_UNMANAGED,

      "Pasture / rangeland" =
        COL_PASTURE
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = x_breaks
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = "Year",
    y = expression(
      "Grassland area (km"^2*")"
    )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.background =
      element_rect(
        fill = "white",
        colour = NA
      ),
    legend.key.size =
      grid::unit(
        0.4,
        "cm"
      ),
    legend.text =
      element_text(
        size = 9
      )
  )

fig <- pA +
  pB +
  plot_annotation(
    tag_levels = "A"
  )

fig_stem <- file.path(
  FIG_DIR,
  "FIG_5_habitat_change_1899_2019"
)

ggsave(
  paste0(fig_stem, ".png"),
  fig,
  width = 11,
  height = 4.2,
  dpi = 320
)

ggsave(
  paste0(fig_stem, ".pdf"),
  fig,
  width = 11,
  height = 4.2
)

ggsave(
  paste0(fig_stem, ".svg"),
  fig,
  width = 11,
  height = 4.2
)

# ============================================================
# 8. REPORT CHECKS
# ============================================================

grass_min <- grass_total |>
  slice_min(
    grassland_km2,
    n = 1,
    with_ties = FALSE
  )

cat(
  "\n==================== HISTORICAL HABITAT CHANGE ====================\n"
)

cat(
  sprintf(
    "Weighted suitable habitat : %.0f -> %.0f km2 (%+.1f%%)\n",
    weighted_1899,
    weighted_2019,
    100 *
      (
        weighted_2019 /
          weighted_1899 -
          1
      )
  )
)

cat(
  sprintf(
    "Total grassland           : %.0f -> %.0f km2 (%+.1f%%)\n",
    grass_1899,
    grass_2019,
    100 *
      (
        grass_2019 /
          grass_1899 -
          1
      )
  )
)

cat(
  sprintf(
    "Grassland minimum         : %.0f km2 in %d\n",
    grass_min$grassland_km2,
    grass_min$year
  )
)

cat(
  sprintf(
    "Terrestrial area 1899     : %.0f km2\n",
    land_totals$land_km2[
      land_totals$year == 1899
    ]
  )
)

cat(
  sprintf(
    "Terrestrial area 2019     : %.0f km2\n",
    land_totals$land_km2[
      land_totals$year == 2019
    ]
  )
)

cat(
  "\n==================== TABLE S5 ====================\n"
)

print(
  as.data.frame(
    table_formatted
  ),
  row.names = FALSE,
  right = FALSE
)

cat(
  "\n==================== OUTPUTS ====================\n"
)

cat(
  "Figure 5 -> ",
  FIG_DIR,
  "/FIG_5_habitat_change_1899_2019.*\n",
  sep = ""
)

cat(
  "Table S5 -> ",
  TAB_DIR,
  "/TABLE_S5_land_cover_change.tsv\n",
  sep = ""
)

cat(
  "Summary -> ",
  TAB_DIR,
  "/habitat_change_summary.tsv\n",
  sep = ""
)
