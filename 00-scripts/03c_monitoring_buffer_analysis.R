# ============================================================
# Monitoring-site landscape composition analysis
#
# Inputs:
#   02-landscapes/derived/sites_with_classes.tsv
#   02-landscapes/derived/site_buffer_composition_long.tsv
#
# Outputs:
#   02-landscapes/derived/site_buffer_composition_wide.tsv
#   02-landscapes/tables/TABLE_S3_landscape_composition_2km.tsv
#   02-landscapes/tables/TABLE_S3_landscape_composition_2km_numeric.tsv
#   02-landscapes/tables/grassland_association_by_radius.tsv
#   02-landscapes/figures/grassland_association_by_radius.*
#
# Usage:
#   Rscript 00-scripts/03c_monitoring_buffer_analysis.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

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

DER <- file.path(REPO_ROOT, "02-landscapes", "derived")
TAB <- file.path(REPO_ROOT, "02-landscapes", "tables")
FIG <- file.path(REPO_ROOT, "02-landscapes", "figures")

dir.create(DER, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

FOCAL_RADIUS    <- 2000
MIN_TERRESTRIAL <- 0.50

GROUP_ORDER <- c(
  "Grassland",
  "Cropland",
  "Forest",
  "Urban",
  "Other"
)

COMP_TERMS <- c(
  "p_grassland",
  "p_cropland",
  "p_forest",
  "p_urban"
)

# ============================================================
# 2. INPUT
# ============================================================

sites <- read_tsv(
  file.path(DER, "sites_with_classes.tsv"),
  show_col_types = FALSE
)

comp <- read_tsv(
  file.path(DER, "site_buffer_composition_long.tsv"),
  show_col_types = FALSE
)

required_site_cols <- c(
  "sit_uid",
  "presence",
  "lat",
  "lon",
  "in_study_region"
)

missing_cols <- setdiff(required_site_cols, names(sites))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in sites_with_classes.tsv: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (any(duplicated(sites$sit_uid))) {
  stop("Duplicated sit_uid values found in sites_with_classes.tsv.")
}

cat(
  "Monitoring sites:", nrow(sites),
  "| recorded presences:", sum(sites$presence, na.rm = TRUE), "\n"
)

cat(
  "Study region:", sum(sites$in_study_region, na.rm = TRUE),
  "| recorded presences:",
  sum(sites$presence[sites$in_study_region], na.rm = TRUE),
  "\n\n"
)

# ============================================================
# 3. CORINE HABITAT GROUPING
# ============================================================

# The buffer analysis uses CORINE land cover only.
# Marine cells are excluded from the terrestrial denominator.

clc_to_group <- function(code) {

  code <- as.integer(code)

  case_when(
    is.na(code) ~ NA_character_,

    # Estuaries / marine water
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

comp_corine <- comp |>
  filter(source == "corine") |>
  mutate(group = clc_to_group(clc_code))

# ============================================================
# 4. BUFFER COMPOSITION
# ============================================================

# Total CORINE cells within each nominal buffer, including marine cells.
buffer_totals <- comp_corine |>
  group_by(sit_uid, radius_m) |>
  summarise(
    buffer_cells = sum(cells),
    .groups = "drop"
  )

# Habitat proportions are calculated relative to terrestrial area.
comp_wide <- comp_corine |>
  filter(!is.na(group)) |>
  group_by(sit_uid, radius_m, group) |>
  summarise(
    cells = sum(cells),
    .groups = "drop"
  ) |>
  group_by(sit_uid, radius_m) |>
  mutate(
    land_cells = sum(cells),
    prop_land = cells / land_cells
  ) |>
  ungroup() |>
  select(
    sit_uid,
    radius_m,
    group,
    prop_land,
    land_cells
  ) |>
  pivot_wider(
    names_from = group,
    values_from = prop_land,
    values_fill = 0
  ) |>
  left_join(
    buffer_totals,
    by = c("sit_uid", "radius_m")
  ) |>
  mutate(
    terrestrial_frac = land_cells / buffer_cells
  )

# Ensure every broad habitat group exists as a column.
for (g in GROUP_ORDER) {
  if (!(g %in% names(comp_wide))) {
    comp_wide[[g]] <- 0
  }
}

comp_wide <- comp_wide |>
  rename_with(
    \(x) paste0("p_", tolower(x)),
    all_of(GROUP_ORDER)
  )

cat(
  "--- buffers excluded for <",
  100 * MIN_TERRESTRIAL,
  "% terrestrial area ---\n",
  sep = ""
)

print(
  as.data.frame(
    comp_wide |>
      group_by(radius_m) |>
      summarise(
        n_total = n(),
        n_excluded = sum(terrestrial_frac < MIN_TERRESTRIAL),
        .groups = "drop"
      )
  )
)

comp_wide <- comp_wide |>
  filter(terrestrial_frac >= MIN_TERRESTRIAL)

# Habitat proportions should sum to one.
stopifnot(
  all(
    abs(
      with(
        comp_wide,
        p_grassland +
          p_cropland +
          p_forest +
          p_urban +
          p_other
      ) - 1
    ) < 1e-8
  )
)

write_tsv(
  comp_wide,
  file.path(DER, "site_buffer_composition_wide.tsv")
)

model_dat <- comp_wide |>
  inner_join(
    sites |>
      select(
        sit_uid,
        presence,
        lat,
        lon,
        in_study_region
      ),
    by = "sit_uid"
  ) |>
  filter(
    !is.na(presence),
    !is.na(lat)
  )

# ============================================================
# 5. TABLE S3 -- LANDSCAPE COMPOSITION WITHIN 2 KM
# ============================================================

dat_reg <- model_dat |>
  filter(
    radius_m == FOCAL_RADIUS,
    in_study_region
  ) |>
  mutate(
    northing_100km =
      (lat - mean(lat, na.rm = TRUE)) / 100000
  )

cat(
  "\n==================== TABLE S3 DATA ====================\n"
)

cat(
  "n =", nrow(dat_reg),
  "| recorded presences =", sum(dat_reg$presence),
  "\n"
)

cat("\nMean habitat cover within 2 km:\n")

print(
  round(
    colMeans(
      dat_reg[, c(COMP_TERMS, "p_other")]
    ),
    3
  )
)

# Other habitat is omitted because habitat proportions sum to one.
m0 <- glm(
  presence ~ northing_100km,
  data = dat_reg,
  family = binomial
)

m1 <- glm(
  presence ~
    p_grassland +
    p_cropland +
    p_forest +
    p_urban +
    northing_100km,
  data = dat_reg,
  family = binomial
)

lrt <- anova(m0, m1, test = "LRT")

or_table <- function(m) {

  s <- summary(m)$coefficients

  tibble(
    term = rownames(s),
    estimate = s[, 1],
    se = s[, 2],
    p = s[, 4]
  ) |>
    mutate(
      scale = if_else(
        term %in% COMP_TERMS,
        0.1,
        1
      ),
      odds_ratio = exp(scale * estimate),
      or_low = exp(scale * (estimate - 1.96 * se)),
      or_high = exp(scale * (estimate + 1.96 * se))
    ) |>
    select(
      term,
      odds_ratio,
      or_low,
      or_high,
      p
    )
}

TERM_LABELS <- c(
  p_grassland = "Grassland cover",
  p_cropland = "Cropland cover",
  p_forest = "Forest cover",
  p_urban = "Urban cover",
  northing_100km = "Northing, per 100 km"
)

table_s3_numeric <- or_table(m1) |>
  filter(term != "(Intercept)") |>
  mutate(
    Predictor = unname(TERM_LABELS[term])
  ) |>
  select(
    Predictor,
    odds_ratio,
    or_low,
    or_high,
    p
  )

write_tsv(
  table_s3_numeric,
  file.path(
    TAB,
    "TABLE_S3_landscape_composition_2km_numeric.tsv"
  )
)

format_p <- function(p) {
  ifelse(
    p < 0.001,
    "<0.001",
    sprintf("%.3f", p)
  )
}

table_s3 <- table_s3_numeric |>
  transmute(
    Predictor,
    `OR per +10 percentage points` =
      sprintf("%.2f", odds_ratio),
    `95% CI` =
      sprintf("%.2f-%.2f", or_low, or_high),
    p = format_p(p)
  )

write_tsv(
  table_s3,
  file.path(
    TAB,
    "TABLE_S3_landscape_composition_2km.tsv"
  )
)

cat(
  "\nLandscape composition beyond northing:",
  "LRT chi2 =", round(lrt$Deviance[2], 2),
  "| df =", lrt$Df[2],
  "| p =", signif(lrt$`Pr(>Chi)`[2], 4),
  "\n\n"
)

print(
  as.data.frame(table_s3),
  right = FALSE
)

# ============================================================
# 6. GRASSLAND ASSOCIATION ACROSS BUFFER RADII
# ============================================================

grassland_by_radius <- function(d) {

  bind_rows(
    lapply(
      sort(unique(d$radius_m)),
      function(rad) {

        dd <- d |>
          filter(radius_m == rad) |>
          mutate(
            northing_100km =
              (lat - mean(lat, na.rm = TRUE)) / 100000
          )

        m <- glm(
          presence ~ p_grassland + northing_100km,
          data = dd,
          family = binomial
        )

        s <- summary(m)$coefficients["p_grassland", ]
        sd_grassland <- sd(dd$p_grassland)

        tibble(
          radius_m = rad,
          n = nrow(dd),
          n_presence = sum(dd$presence),
          mean_grassland = mean(dd$p_grassland),
          sd_grassland = sd_grassland,

          or_per_10pct = exp(0.1 * s[1]),
          or10_low =
            exp(0.1 * (s[1] - 1.96 * s[2])),
          or10_high =
            exp(0.1 * (s[1] + 1.96 * s[2])),

          or_per_sd =
            exp(sd_grassland * s[1]),
          orsd_low =
            exp(
              sd_grassland *
                (s[1] - 1.96 * s[2])
            ),
          orsd_high =
            exp(
              sd_grassland *
                (s[1] + 1.96 * s[2])
            ),

          p = s[4]
        )
      }
    )
  )
}

radius_tab <- grassland_by_radius(
  model_dat |>
    filter(in_study_region)
)

write_tsv(
  radius_tab,
  file.path(
    TAB,
    "grassland_association_by_radius.tsv"
  )
)

cat(
  "\n--- grassland association across buffer radii ---\n"
)

print(
  as.data.frame(
    radius_tab |>
      mutate(
        across(
          where(is.numeric),
          \(x) signif(x, 3)
        )
      )
  )
)

# Supporting plot. This is not numbered as a supplementary figure.
p_radius <- ggplot(
  radius_tab,
  aes(
    x = or_per_sd,
    y = factor(radius_m / 1000)
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    colour = "grey55"
  ) +
  geom_errorbar(
    aes(
      xmin = orsd_low,
      xmax = orsd_high
    ),
    orientation = "y",
    width = 0,
    colour = "grey40"
  ) +
  geom_point(
    size = 2.6,
    colour = "#1B5E20"
  ) +
  scale_x_log10() +
  labs(
    x = "Odds ratio per 1 SD increase in grassland cover",
    y = "Buffer radius (km)"
  ) +
  theme_classic(base_size = 11)

ggsave(
  file.path(
    FIG,
    "grassland_association_by_radius.png"
  ),
  p_radius,
  width = 6.5,
  height = 3.6,
  dpi = 320
)

ggsave(
  file.path(
    FIG,
    "grassland_association_by_radius.pdf"
  ),
  p_radius,
  width = 6.5,
  height = 3.6
)

ggsave(
  file.path(
    FIG,
    "grassland_association_by_radius.svg"
  ),
  p_radius,
  width = 6.5,
  height = 3.6
)

# ============================================================
# 7. OUTPUT SUMMARY
# ============================================================

cat("\n==================== OUTPUTS ====================\n")

cat(
  "Table S3:\n",
  "  TABLE_S3_landscape_composition_2km.tsv\n",
  "  TABLE_S3_landscape_composition_2km_numeric.tsv\n",
  sep = ""
)

cat(
  "\nBuffer-radius sensitivity:\n",
  "  grassland_association_by_radius.tsv\n",
  "  grassland_association_by_radius.*\n",
  sep = ""
)

cat(
  "\nDerived data:\n",
  "  site_buffer_composition_wide.tsv\n",
  sep = ""
)
