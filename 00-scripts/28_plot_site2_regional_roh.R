#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
})

# ============================================================
# PATHS
# ============================================================

BASE <- path.expand(
    "~/tabea_work/08_population_analysis_final/skane_only/results/regional_roh_site2"
)

COMBINED <- file.path(
    BASE,
    "combined"
)

OUT <- file.path(
    BASE,
    "plots"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# INPUT
# ============================================================

length_file <- file.path(
    COMBINED,
    "regional_ROH_length_class_by_replicate.tsv"
)

individual_file <- file.path(
    COMBINED,
    "regional_ROH_all_individuals.tsv"
)

length_rep <- read_tsv(
    length_file,
    show_col_types = FALSE
)

individual <- read_tsv(
    individual_file,
    show_col_types = FALSE
)

# ============================================================
# ORDERS / LABELS
# ============================================================

state_order <- c(
    "1900",
    "2020",
    "2140 SQ",
    "2140 R2",
    "2140 R4",
    "2140 R6"
)

state_labels <- c(
    "1900",
    "2020",
    "2140\nStatus quo",
    "2140\nRestore 2 km",
    "2140\nRestore 4 km",
    "2140\nRestore 6 km"
)

region_order <- c(
    "NW",
    "NE",
    "W",
    "E",
    "SE"
)

roh_order <- c(
    "100-250kb",
    "250-500kb",
    "500kb-1Mb",
    ">1Mb"
)

roh_labels <- c(
    "100–250 kb",
    "250–500 kb",
    "500 kb–1 Mb",
    ">1 Mb"
)

# Darker = longer ROH
roh_cols <- c(
    "100-250kb" = "#D9E6F2",
    "250-500kb" = "#91B8D8",
    "500kb-1Mb" = "#4F83B6",
    ">1Mb"       = "#194F7A"
)

# ============================================================
# COMMON THEME
# ============================================================

theme_roh <- theme_bw(base_size = 11) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),

        strip.background = element_rect(
            fill = "white",
            colour = "black",
            linewidth = 0.4
        ),

        strip.text = element_text(
            face = "bold",
            size = 10
        ),

        axis.title = element_text(
            size = 11
        ),

        axis.text = element_text(
            size = 9
        ),

        legend.title = element_blank(),

        legend.position = "bottom",

        legend.text = element_text(
            size = 9
        ),

        plot.title = element_text(
            face = "bold",
            size = 12
        ),

        plot.margin = margin(
            8, 8, 8, 8
        )
    )

# ============================================================
# FIGURE 1
#
# POPULATION-LEVEL ROH COMPOSITION
#
# One panel = one state x region.
# Each stacked bar = mean FROH contribution across replicates.
# ============================================================

pop_summary <- length_rep %>%
    mutate(
        state = factor(
            state,
            levels = state_order,
            labels = state_labels
        ),

        region = factor(
            region,
            levels = region_order
        ),

        roh_class = factor(
            roh_class,
            levels = roh_order,
            labels = roh_labels
        )
    ) %>%
    group_by(
        state,
        region,
        roh_class
    ) %>%
    summarise(
        mean_FROH = mean(
            FROH_contribution,
            na.rm = TRUE
        ),

        sd_FROH = sd(
            FROH_contribution,
            na.rm = TRUE
        ),

        n_reps = sum(
            is.finite(FROH_contribution)
        ),

        .groups = "drop"
    )

# Save exact plotted values
write_tsv(
    pop_summary,
    file.path(
        OUT,
        "Figure_ROH_population_composition_values.tsv"
    )
)

p_population <- ggplot(
    pop_summary,
    aes(
        x = 1,
        y = mean_FROH,
        fill = roh_class
    )
) +

    geom_col(
        width = 0.68
    ) +

    facet_grid(
        state ~ region,
        switch = "y"
    ) +

    scale_fill_manual(
        values = setNames(
            roh_cols,
            roh_labels
        ),
        drop = FALSE
    ) +

    scale_x_continuous(
        breaks = NULL
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0, 0.05)
        )
    ) +

    labs(
        x = NULL,
        y = expression(F[ROH]),
        fill = NULL
    ) +

    theme_roh +

    theme(
        strip.placement = "outside",

        strip.text.y.left = element_text(
            angle = 0
        ),

        panel.spacing.x = unit(
            0.7,
            "lines"
        ),

        panel.spacing.y = unit(
            0.35,
            "lines"
        )
    )

ggsave(
    file.path(
        OUT,
        "Figure_ROH_population_composition.svg"
    ),
    p_population,
    width = 9.0,
    height = 9.0,
    units = "in"
)

ggsave(
    file.path(
        OUT,
        "Figure_ROH_population_composition.png"
    ),
    p_population,
    width = 9.0,
    height = 9.0,
    units = "in",
    dpi = 400
)

# ============================================================
# FIGURE 2
#
# INDIVIDUAL ROH PROFILES
#
# One stacked bar = one individual.
# Individuals are ordered within state x region by total FROH.
#
# This shows whether regional patterns are widespread among
# individuals or driven by a small number of highly inbred
# individuals.
# ============================================================

individual_long <- individual %>%
    mutate(
        state = case_when(
            scenario %in% c(
                "historical",
                "shared_history"
            ) &
                year == 1900 ~ "1900",

            scenario %in% c(
                "historical",
                "shared_history"
            ) &
                year == 2020 ~ "2020",

            scenario == "status_quo" ~
                "2140 SQ",

            scenario == "restore_2km" ~
                "2140 R2",

            scenario == "restore_4km" ~
                "2140 R4",

            scenario == "restore_6km" ~
                "2140 R6",

            TRUE ~ NA_character_
        )
    ) %>%

    filter(
        !is.na(state)
    ) %>%

    mutate(
        state = factor(
            state,
            levels = state_order,
            labels = state_labels
        ),

        region = factor(
            region,
            levels = region_order
        )
    ) %>%

    group_by(
        state,
        region
    ) %>%

    arrange(
        FROH_100kb,
        .by_group = TRUE
    ) %>%

    mutate(
        individual_order = row_number()
    ) %>%

    ungroup() %>%

    select(
        state,
        region,
        rep,
        site,
        ind_index,
        individual_order,
        FROH_100kb,
        FROH_100_250kb,
        FROH_250_500kb,
        FROH_500kb_1Mb,
        FROH_gt1Mb
    ) %>%

    pivot_longer(
        cols = c(
            FROH_100_250kb,
            FROH_250_500kb,
            FROH_500kb_1Mb,
            FROH_gt1Mb
        ),

        names_to = "roh_class",
        values_to = "FROH_contribution"
    ) %>%

    mutate(
        roh_class = recode(
            roh_class,
            "FROH_100_250kb" =
                "100–250 kb",

            "FROH_250_500kb" =
                "250–500 kb",

            "FROH_500kb_1Mb" =
                "500 kb–1 Mb",

            "FROH_gt1Mb" =
                ">1 Mb"
        ),

        roh_class = factor(
            roh_class,
            levels = roh_labels
        )
    )

write_tsv(
    individual_long,
    file.path(
        OUT,
        "Figure_ROH_individual_profiles_values.tsv"
    )
)

p_individual <- ggplot(
    individual_long,
    aes(
        x = individual_order,
        y = FROH_contribution,
        fill = roh_class
    )
) +

    geom_col(
        width = 1
    ) +

    facet_grid(
        state ~ region,
        scales = "free_x",
        switch = "y"
    ) +

    scale_fill_manual(
        values = setNames(
            roh_cols,
            roh_labels
        ),
        drop = FALSE
    ) +

    scale_x_continuous(
        breaks = NULL,
        expand = expansion(
            mult = c(0, 0)
        )
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0, 0.03)
        )
    ) +

    labs(
        x = "Individuals",
        y = expression(F[ROH]),
        fill = NULL
    ) +

    theme_roh +

    theme(
        strip.placement = "outside",

        strip.text.y.left = element_text(
            angle = 0
        ),

        panel.spacing.x = unit(
            0.35,
            "lines"
        ),

        panel.spacing.y = unit(
            0.35,
            "lines"
        )
    )

ggsave(
    file.path(
        OUT,
        "Figure_ROH_individual_profiles.svg"
    ),
    p_individual,
    width = 12,
    height = 9,
    units = "in"
)

ggsave(
    file.path(
        OUT,
        "Figure_ROH_individual_profiles.png"
    ),
    p_individual,
    width = 12,
    height = 9,
    units = "in",
    dpi = 400
)

# ============================================================
# OPTIONAL: ONE INDIVIDUAL-PROFILE FIGURE PER REGION
#
# Much easier to read than all five regions together.
# ============================================================

for (reg in region_order) {

    dat_reg <- individual_long %>%
        filter(
            region == reg
        )

    p_reg <- ggplot(
        dat_reg,
        aes(
            x = individual_order,
            y = FROH_contribution,
            fill = roh_class
        )
    ) +

        geom_col(
            width = 1
        ) +

        facet_wrap(
            ~ state,
            ncol = 1,
            scales = "free_x"
        ) +

        scale_fill_manual(
            values = setNames(
                roh_cols,
                roh_labels
            ),
            drop = FALSE
        ) +

        scale_x_continuous(
            breaks = NULL,
            expand = expansion(
                mult = c(0, 0)
            )
        ) +

        scale_y_continuous(
            expand = expansion(
                mult = c(0, 0.03)
            )
        ) +

        labs(
            title = paste(
                "Region",
                reg
            ),
            x = "Individuals",
            y = expression(F[ROH]),
            fill = NULL
        ) +

        theme_roh

    ggsave(
        file.path(
            OUT,
            paste0(
                "Figure_ROH_individual_profiles_",
                reg,
                ".svg"
            )
        ),
        p_reg,
        width = 7,
        height = 9,
        units = "in"
    )

    ggsave(
        file.path(
            OUT,
            paste0(
                "Figure_ROH_individual_profiles_",
                reg,
                ".png"
            )
        ),
        p_reg,
        width = 7,
        height = 9,
        units = "in",
        dpi = 400
    )
}

# ============================================================
# FINISH
# ============================================================

cat("\n============================================\n")
cat("ROH FIGURES COMPLETE\n")
cat("============================================\n\n")

cat("Population-level figure:\n")
cat(
    file.path(
        OUT,
        "Figure_ROH_population_composition.svg"
    ),
    "\n\n"
)

cat("Individual-level figure:\n")
cat(
    file.path(
        OUT,
        "Figure_ROH_individual_profiles.svg"
    ),
    "\n\n"
)

cat("Individual regional figures:\n")

for (reg in region_order) {
    cat(
        "  ",
        file.path(
            OUT,
            paste0(
                "Figure_ROH_individual_profiles_",
                reg,
                ".svg"
            )
        ),
        "\n"
    )
}

cat("\nDONE\n")
