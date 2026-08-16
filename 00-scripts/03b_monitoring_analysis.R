# ============================================================
# Monitoring-site habitat association analysis
#
# Inputs:
#   02-landscapes/derived/sites_with_classes.tsv
#
# Outputs:
#   02-landscapes/tables/TABLE_S1_habitat_occurrence.tsv
#   02-landscapes/tables/TABLE_S1_habitat_occurrence_numeric.tsv
#   02-landscapes/tables/TABLE_S2_logistic_regression.tsv
#   02-landscapes/tables/TABLE_S2_logistic_regression_numeric.tsv
#   02-landscapes/tables/FIG_S2_CORINE_HILDA_correspondence_data.tsv
#   02-landscapes/figures/FIG_S1_habitat_association.*
#   02-landscapes/figures/FIG_S2_CORINE_HILDA_correspondence.*
#
# Usage:
#   Rscript 00-scripts/03b_monitoring_analysis.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

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

dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

MIN_SITES <- 10

# ==============================================================================
# 1. DATA AND HABITAT RECLASSIFICATION
# ==============================================================================

sites <- read_tsv(
  file.path(DER, "sites_with_classes.tsv"),
  show_col_types = FALSE
)

required_cols <- c(
  "presence",
  "in_study_region",
  "lat",
  "clc_code_extracted",
  "clc_name_extracted",
  "hilda_class"
)

missing_cols <- setdiff(required_cols, names(sites))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in sites_with_classes.tsv: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Explicit CORINE -> broad habitat lookup.
# Sparse/bare land, wetlands and water are pooled into Other.
corine_lookup <- tibble(
  clc_code = c(
    231,
    321, 322, 324,
    211, 242, 243,
    311, 312, 313,
    112, 121, 122, 124, 131, 141, 142,
    331, 332, 333, 334, 411, 412, 512, 523
  ),
  corine_class = c(
    "Pasture/Rangeland",
    rep("Unmanaged grass/shrubland", 3),
    rep("Cropland", 3),
    rep("Forest", 3),
    rep("Urban", 7),
    rep("Other", 8)
  )
)

CLASS_ORDER <- c(
  "Pasture/Rangeland",
  "Unmanaged grass/shrubland",
  "Cropland",
  "Forest",
  "Urban",
  "Other"
)

ANALYSIS_ORDER <- c(
  "Grassland",
  "Cropland",
  "Forest",
  "Urban",
  "Other"
)

# Use the CORINE class extracted directly from the 2018 raster in 03a.
# The extraction reproduced the existing site classifications with 100% agreement.
sites <- sites |>
  mutate(
    clc_code = as.integer(clc_code_extracted),
    clc_name = clc_name_extracted
  ) |>
  left_join(corine_lookup, by = "clc_code")

# Stop if a CORINE class has not been explicitly mapped.
unmapped <- sites |>
  filter(is.na(corine_class)) |>
  distinct(clc_code, clc_name)

if (nrow(unmapped) > 0) {
  print(unmapped)
  stop("Unmapped CORINE classes found. Add them to corine_lookup before continuing.")
}

# Six descriptive classes + pooled Grassland class for regression/figures.
sites <- sites |>
  mutate(
    corine_class = factor(corine_class, levels = CLASS_ORDER),
    analysis_class = if_else(
      as.character(corine_class) %in%
        c("Pasture/Rangeland", "Unmanaged grass/shrubland"),
      "Grassland",
      as.character(corine_class)
    ),
    analysis_class = factor(analysis_class, levels = ANALYSIS_ORDER)
  )

# Reclassify HILDA+ to the same broad categories.
hilda_to_broad <- function(x) {
  z <- tolower(trimws(as.character(x)))

  case_when(
    is.na(x)                               ~ NA_character_,
    z == "pasture / rangeland"            ~ "Pasture/Rangeland",
    z == "unmanaged grass / shrubland"    ~ "Unmanaged grass/shrubland",
    z == "cropland"                       ~ "Cropland",
    z == "forest"                         ~ "Forest",
    z == "urban"                          ~ "Urban",
    TRUE                                   ~ "Other"
  )
}

sites <- sites |>
  mutate(
    hilda_broad_class = hilda_to_broad(hilda_class),
    hilda_analysis_class = if_else(
      hilda_broad_class %in%
        c("Pasture/Rangeland", "Unmanaged grass/shrubland"),
      "Grassland",
      hilda_broad_class
    )
  )

# Overall recorded-presence proportions used for Figure S1.
overall_all <- mean(sites$presence, na.rm = TRUE)
overall_reg <- mean(
  sites$presence[sites$in_study_region],
  na.rm = TRUE
)

cat("\n==================== DATA SUMMARY ====================\n")
cat(
  "Broader monitoring area:", nrow(sites),
  "| sites with recorded presence:", sum(sites$presence, na.rm = TRUE),
  "| proportion:", round(overall_all, 3), "\n"
)
cat(
  "Study region:", sum(sites$in_study_region, na.rm = TRUE),
  "| sites with recorded presence:",
  sum(sites$presence[sites$in_study_region], na.rm = TRUE),
  "| proportion:", round(overall_reg, 3), "\n"
)

# CORINE codes collapsed into one string for Table S1.
class_codes <- corine_lookup |>
  mutate(corine_class = factor(corine_class, levels = CLASS_ORDER)) |>
  group_by(corine_class) |>
  summarise(
    clc_codes = paste(sort(clc_code), collapse = ", "),
    .groups = "drop"
  ) |>
  rename(class = corine_class)

# ==============================================================================
# 2. HELPERS
# ==============================================================================

# Proportion of monitoring sites with recorded presence.
# Confidence intervals are exact Clopper-Pearson 95% binomial intervals.
proportion_summary <- function(d, class_col) {
  d |>
    filter(!is.na(.data[[class_col]]), !is.na(presence)) |>
    group_by(class = .data[[class_col]]) |>
    summarise(
      n_sites = n(),
      n_recorded = sum(presence),
      .groups = "drop"
    ) |>
    rowwise() |>
    mutate(
      proportion = n_recorded / n_sites,
      ci_low = binom.test(n_recorded, n_sites)$conf.int[1],
      ci_high = binom.test(n_recorded, n_sites)$conf.int[2]
    ) |>
    ungroup()
}

# Logistic-regression coefficient table.
# Confidence intervals are Wald 95% intervals on the log-odds scale,
# exponentiated to give odds-ratio intervals.
or_table <- function(m) {
  s <- summary(m)$coefficients

  data.frame(
    term = rownames(s),
    estimate = s[, 1],
    se = s[, 2],
    p = s[, 4],
    row.names = NULL
  ) |>
    mutate(
      odds_ratio = exp(estimate),
      or_low = exp(estimate - 1.96 * se),
      or_high = exp(estimate + 1.96 * se)
    ) |>
    select(term, odds_ratio, or_low, or_high, p)
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "-",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

# ==============================================================================
# 3. DESCRIPTIVE HABITAT SUMMARIES
# ==============================================================================

# Six classes are kept separate for Table S1.
broad_all <- proportion_summary(sites, "corine_class")
broad_reg <- proportion_summary(
  sites |> filter(in_study_region),
  "corine_class"
)

# Grassland is pooled for Figure S1 and the regression.
primary_all <- proportion_summary(sites, "analysis_class")
primary_reg <- proportion_summary(
  sites |> filter(in_study_region),
  "analysis_class"
)

# ==============================================================================
# 4. TABLE S1 -- HABITAT CLASSIFICATION AND MONITORING-SITE OCCURRENCE
# ==============================================================================

side_desc <- function(x, suffix) {
  x |>
    mutate(class = as.character(class)) |>
    select(class, n_sites, n_recorded, proportion, ci_low, ci_high) |>
    rename_with(\(z) paste0(z, suffix), -class)
}

table_s1_numeric <- class_codes |>
  mutate(class = as.character(class)) |>
  left_join(side_desc(broad_reg, "_reg"), by = "class") |>
  left_join(side_desc(broad_all, "_all"), by = "class") |>
  mutate(class = factor(class, levels = CLASS_ORDER)) |>
  arrange(class)

write_tsv(
  table_s1_numeric,
  file.path(TAB, "TABLE_S1_habitat_occurrence_numeric.tsv")
)

fmt_prop <- function(p, lo, hi) {
  ifelse(
    is.na(p),
    "-",
    sprintf("%.1f (%.1f-%.1f)", 100 * p, 100 * lo, 100 * hi)
  )
}

table_s1 <- table_s1_numeric |>
  transmute(
    `Habitat class` = class,
    `CORINE codes` = clc_codes,
    `Monitoring sites (study region)` = n_sites_reg,
    `Sites with recorded presence (study region)` = n_recorded_reg,
    `Recorded presence, % (95% CI), study region` =
      fmt_prop(proportion_reg, ci_low_reg, ci_high_reg),
    `Monitoring sites (broader area)` = n_sites_all,
    `Sites with recorded presence (broader area)` = n_recorded_all,
    `Recorded presence, % (95% CI), broader area` =
      fmt_prop(proportion_all, ci_low_all, ci_high_all)
  )

write_tsv(
  table_s1,
  file.path(TAB, "TABLE_S1_habitat_occurrence.tsv")
)

cat("\n==================== TABLE S1 ====================\n")
print(as.data.frame(table_s1), right = FALSE)

# ==============================================================================
# 5. TABLE S2 -- LOGISTIC REGRESSION OF RECORDED PRESENCE
# ==============================================================================

run_habitat_model <- function(d, label) {

  # Only groups with >= MIN_SITES enter the model.
  keep <- d |>
    filter(!is.na(analysis_class)) |>
    count(analysis_class) |>
    filter(n >= MIN_SITES) |>
    pull(analysis_class)

  dd <- d |>
    filter(
      analysis_class %in% keep,
      !is.na(presence),
      !is.na(lat)
    ) |>
    mutate(
      analysis_ref = relevel(
        droplevels(factor(analysis_class)),
        ref = "Grassland"
      ),
      # 'lat' in the input is SWEREF99 TM northing.
      # Coefficient is therefore interpreted per 100 km northward.
      northing_100km = (lat - mean(lat, na.rm = TRUE)) / 100000
    )

  if (!("Grassland" %in% levels(dd$analysis_ref))) {
    stop("Grassland is unavailable as the reference class in ", label)
  }

  # Null and habitat models use exactly the same observations.
  m0 <- glm(
    presence ~ northing_100km,
    data = dd,
    family = binomial
  )

  m1 <- glm(
    presence ~ analysis_ref + northing_100km,
    data = dd,
    family = binomial
  )

  lrt <- anova(m0, m1, test = "LRT")
  coefs <- or_table(m1)

  habitat_terms <- coefs |>
    filter(grepl("^analysis_ref", term)) |>
    mutate(predictor = sub("^analysis_ref", "", term)) |>
    select(predictor, odds_ratio, or_low, or_high, p)

  habitat_terms <- tibble(predictor = ANALYSIS_ORDER) |>
    left_join(habitat_terms, by = "predictor") |>
    mutate(
      odds_ratio = if_else(predictor == "Grassland", 1, odds_ratio),
      or_low = if_else(predictor == "Grassland", NA_real_, or_low),
      or_high = if_else(predictor == "Grassland", NA_real_, or_high),
      p = if_else(predictor == "Grassland", NA_real_, p)
    )

  northing_term <- coefs |>
    filter(term == "northing_100km") |>
    transmute(
      predictor = "Northing (per 100 km)",
      odds_ratio,
      or_low,
      or_high,
      p
    )

  list(
    coefficients = bind_rows(habitat_terms, northing_term),
    overall = tibble(
      extent = label,
      n_modelled = nrow(dd),
      lrt_chisq = unname(lrt$Deviance[2]),
      lrt_df = unname(lrt$Df[2]),
      lrt_p = unname(lrt$`Pr(>Chi)`[2])
    )
  )
}

model_reg <- run_habitat_model(
  sites |> filter(in_study_region),
  "Study region"
)

model_all <- run_habitat_model(
  sites,
  "Broader monitoring area"
)

side_or <- function(x, suffix) {
  x |>
    rename_with(\(z) paste0(z, suffix), -predictor)
}

predictor_order <- c(
  ANALYSIS_ORDER,
  "Northing (per 100 km)"
)

table_s2_numeric <- tibble(predictor = predictor_order) |>
  left_join(side_or(model_reg$coefficients, "_reg"), by = "predictor") |>
  left_join(side_or(model_all$coefficients, "_all"), by = "predictor")

write_tsv(
  table_s2_numeric,
  file.path(TAB, "TABLE_S2_logistic_regression_numeric.tsv")
)

fmt_or <- function(predictor, est, lo, hi) {
  ifelse(
    predictor == "Grassland",
    "1.00 (reference)",
    sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
  )
}

table_s2 <- table_s2_numeric |>
  transmute(
    Predictor = predictor,
    `Study region OR (95% CI)` =
      fmt_or(predictor, odds_ratio_reg, or_low_reg, or_high_reg),
    `Study region p` = format_p(p_reg),
    `Broader area OR (95% CI)` =
      fmt_or(predictor, odds_ratio_all, or_low_all, or_high_all),
    `Broader area p` = format_p(p_all)
  )

write_tsv(
  table_s2,
  file.path(TAB, "TABLE_S2_logistic_regression.tsv")
)

# Overall habitat-effect tests are used in the Table S2 caption, not as a
# separate supplementary table.
overall_tests <- bind_rows(model_reg$overall, model_all$overall)

cat("\n==================== TABLE S2 ====================\n")
print(as.data.frame(table_s2), right = FALSE)

cat("\nOverall habitat effect (for Table S2 caption):\n")
print(
  overall_tests |>
    mutate(
      lrt_chisq = round(lrt_chisq, 2),
      lrt_p = signif(lrt_p, 4)
    )
)

# ==============================================================================
# 6. FIGURE S1 -- RECORDED PRESENCE BY BROAD HABITAT GROUP
# ==============================================================================

primary_both <- bind_rows(
  primary_reg |> mutate(extent = "Study region"),
  primary_all |> mutate(extent = "Broader monitoring area")
) |>
  mutate(
    class = factor(as.character(class), levels = rev(ANALYSIS_ORDER)),
    extent = factor(
      extent,
      levels = c("Study region", "Broader monitoring area")
    )
  )

overall_lines <- tibble(
  extent = factor(
    c("Study region", "Broader monitoring area"),
    levels = c("Study region", "Broader monitoring area")
  ),
  value = c(overall_reg, overall_all)
)

p_habitat <- ggplot(
  primary_both,
  aes(x = proportion, y = class)
) +
  geom_vline(
    data = overall_lines,
    aes(xintercept = value),
    linetype = 2,
    colour = "grey55",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0,
    colour = "grey45",
    linewidth = 0.7
  ) +
  geom_point(
    colour = "#1B5E20",
    size = 2.8
  ) +
  geom_text(
    aes(
      x = ci_high + 0.012,
      label = paste0("n = ", n_sites)
    ),
    hjust = 0,
    size = 3.1,
    colour = "grey35"
  ) +
  facet_wrap(~extent, nrow = 1) +
  scale_x_continuous(
    limits = c(0, 0.55),
    breaks = seq(0, 0.5, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = expression(
      paste(
        "Proportion of monitoring sites with ",
        italic("C. semiargus"),
        " recorded"
      )
    ),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.spacing = grid::unit(1.6, "lines"),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
  )

ggsave(
  file.path(FIG, "FIG_S1_habitat_association.png"),
  p_habitat,
  width = 9,
  height = 4.5,
  dpi = 320
)

ggsave(
  file.path(FIG, "FIG_S1_habitat_association.pdf"),
  p_habitat,
  width = 9,
  height = 4.5
)

ggsave(
  file.path(FIG, "FIG_S1_habitat_association.svg"),
  p_habitat,
  width = 9,
  height = 4.5
)

# ==============================================================================
# 7. FIGURE S2 -- CORINE -> HILDA+ HABITAT CORRESPONDENCE
# Study region only; grasslands pooled; each CORINE row sums to 100%.
# ==============================================================================

heatmap_dat <- sites |>
  filter(
    in_study_region,
    !is.na(analysis_class),
    !is.na(hilda_analysis_class)
  ) |>
  mutate(
    CORINE = factor(
      as.character(analysis_class),
      levels = ANALYSIS_ORDER
    ),
    HILDA = factor(
      hilda_analysis_class,
      levels = ANALYSIS_ORDER
    )
  ) |>
  count(CORINE, HILDA, name = "n") |>
  complete(
    CORINE,
    HILDA,
    fill = list(n = 0)
  ) |>
  group_by(CORINE) |>
  mutate(
    row_total = sum(n),
    proportion = n / row_total,
    percent = 100 * proportion,
    label = if_else(
      n == 0,
      "",
      paste0(
        sprintf("%.1f", percent), "%\n",
        "(n = ", n, ")"
      )
    )
  ) |>
  ungroup()

# Underlying Figure S2 values are retained for reproducibility.
write_tsv(
  heatmap_dat,
  file.path(TAB, "FIG_S2_CORINE_HILDA_correspondence_data.tsv")
)

# Sanity check: rows should sum to 100%.
cat("\nFigure S2 row-normalization check:\n")
print(
  heatmap_dat |>
    group_by(CORINE) |>
    summarise(
      n_sites = sum(n),
      percent_total = sum(percent),
      .groups = "drop"
    )
)

p_corine_hilda <- ggplot(
  heatmap_dat,
  aes(x = HILDA, y = CORINE, fill = proportion)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = label),
    size = 3.3,
    lineheight = 0.9
  ) +
  scale_fill_gradient(
    low = "#F3F7F2",
    high = "#1B5E20",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    name = "Within-CORINE\nclass (%)"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(
    drop = FALSE,
    limits = rev(ANALYSIS_ORDER)
  ) +
  labs(
    x = "HILDA+ habitat class (~1 km)",
    y = "CORINE habitat class (100 m)"
  ) +
  coord_fixed() +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.title = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(FIG, "FIG_S2_CORINE_HILDA_correspondence.png"),
  p_corine_hilda,
  width = 7.2,
  height = 5.8,
  dpi = 320
)

ggsave(
  file.path(FIG, "FIG_S2_CORINE_HILDA_correspondence.pdf"),
  p_corine_hilda,
  width = 7.2,
  height = 5.8
)

ggsave(
  file.path(FIG, "FIG_S2_CORINE_HILDA_correspondence.svg"),
  p_corine_hilda,
  width = 7.2,
  height = 5.8
)

# ==============================================================================
# 8. OUTPUT SUMMARY
# ==============================================================================

cat("\n==================== OUTPUTS ====================\n")
cat("Table S1:\n")
cat("  TABLE_S1_habitat_occurrence.tsv\n")
cat("  TABLE_S1_habitat_occurrence_numeric.tsv\n")
cat("\nTable S2:\n")
cat("  TABLE_S2_logistic_regression.tsv\n")
cat("  TABLE_S2_logistic_regression_numeric.tsv\n")
cat("  Overall LRT statistics are printed above for the Table S2 caption.\n")
cat("\nFigure S1:\n")
cat("  FIG_S1_habitat_association.*\n")
cat("\nFigure S2:\n")
cat("  FIG_S2_CORINE_HILDA_correspondence.*\n")
cat("  FIG_S2_CORINE_HILDA_correspondence_data.tsv\n")
cat("\nTables ->", TAB, "\nFigures ->", FIG, "\n")

# Display figures when the script is run interactively.
p_habitat
p_corine_hilda