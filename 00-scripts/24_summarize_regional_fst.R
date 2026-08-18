#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

INFILE <- if (length(args) >= 1) {
  args[1]
} else {
  "08-population-analysis/full/results/regional_genetics/regional_fst_by_replicate.tsv"
}

OUTDIR <- if (length(args) >= 2) {
  args[2]
} else {
  "08-population-analysis/full/results/regional_genetics"
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

SKANE <- c("NW", "NE", "W", "E", "SE")

fst <- read_tsv(INFILE, show_col_types = FALSE) %>%

  mutate(
    fst_group = case_when(

      region1 %in% SKANE &
        region2 %in% SKANE ~
        "Within Skane",

      (region1 == "NEW02" & region2 %in% SKANE) |
        (region2 == "NEW02" & region1 %in% SKANE) ~
        "Skane - W Smaland",

      (region1 == "NEW01" & region2 %in% SKANE) |
        (region2 == "NEW01" & region1 %in% SKANE) ~
        "Skane - Oland",

      TRUE ~ NA_character_
    ),

    state = case_when(

      scenario == "shared_history" & year == 1900 ~
        "1900",

      scenario == "shared_history" & year == 2020 ~
        "2020",

      scenario == "status_quo" & year == 2140 ~
        "2140 status quo",

      scenario == "restore_2km" & year == 2140 ~
        "2140 restore 2 km",

      scenario == "restore_4km" & year == 2140 ~
        "2140 restore 4 km",

      scenario == "restore_6km" & year == 2140 ~
        "2140 restore 6 km",

      TRUE ~ NA_character_
    )
  ) %>%

  filter(
    !is.na(fst_group),
    !is.na(state)
  )

# First average pairwise FST values within each simulation replicate.
fst_rep <- fst %>%
  group_by(state, fst_group, rep) %>%
  summarise(
    mean_fst_rep = if (
      all(is.na(fst_hudson))
    ) {
      NA_real_
    } else {
      mean(fst_hudson, na.rm = TRUE)
    },
    n_pairs = sum(!is.na(fst_hudson)),
    .groups = "drop"
  )

# Then summarize across independent simulation replicates.
fst_summary <- fst_rep %>%
  group_by(state, fst_group) %>%
  summarise(
    n_reps = sum(!is.na(mean_fst_rep)),
    mean_fst = mean(mean_fst_rep, na.rm = TRUE),
    sd_fst = sd(mean_fst_rep, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  fst_rep,
  file.path(OUTDIR, "fst_group_means_by_replicate.tsv")
)

write_tsv(
  fst_summary,
  file.path(OUTDIR, "fst_group_summary.tsv")
)

print(fst_summary, n = Inf)
