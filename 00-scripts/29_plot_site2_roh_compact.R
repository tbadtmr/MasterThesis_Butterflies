#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
})

# ============================================================
# PATHS
# ============================================================

BASE <- path.expand(
    "~/tabea_work/08_population_analysis_final/skane_only/results/regional_roh_site2"
)

INFILE <- file.path(
    BASE,
    "combined",
    "regional_ROH_length_class_by_replicate.tsv"
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

dat <- read_tsv(
    INFILE,
    show_col_types = FALSE
)

# ============================================================
# ORDER AND LABELS
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
    "SQ",
    "R2",
    "R4",
    "R6"
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

# ============================================================
# COLOURS
# ============================================================

# ROH classes:
# restrained grey scale similar to the empirical ROH figure.
# Darker grey = longer ROH.

roh_colours <- c(
    "100–250 kb" = "#D9D9D9",
    "250–500 kb" = "#B7B7B7",
    "500 kb–1 Mb" = "#858585",
    ">1 Mb" = "#4A4A4A"
)

# Scenario colours consistent with the other thesis figures.
# These are used only for replicate points and summary markers.

scenario_colours <- c(
    "1900" = "#9A9A9A",
    "2020" = "#666666",
    "SQ"   = "#000000",
    "R2"   = "#69A541",
    "R4"   = "#21A6B6",
    "R6"   = "#9A6DCE"
)

# ============================================================
# PREPARE DATA
# ============================================================

dat <- dat %>%
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
    )

# ============================================================
# MEAN ROH-LENGTH COMPOSITION ACROSS REPLICATES
# ============================================================

class_mean <- dat %>%
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
        .groups = "drop"
    )

# ============================================================
# TOTAL FROH FOR EACH SIMULATION REPLICATE
# ============================================================

rep_total <- dat %>%
    group_by(
        state,
        region,
        rep
    ) %>%
    summarise(
        total_FROH = sum(
            FROH_contribution,
            na.rm = TRUE
        ),
        .groups = "drop"
    )

# ============================================================
# MEAN ± SD ACROSS REPLICATES
# ============================================================

total_summary <- rep_total %>%
    group_by(
        state,
        region
    ) %>%
    summarise(
        mean_total = mean(
            total_FROH,
            na.rm = TRUE
        ),

        sd_total = sd(
            total_FROH,
            na.rm = TRUE
        ),

        n_reps = n(),

        .groups = "drop"
    )

write_tsv(
    total_summary,
    file.path(
        OUT,
        "Figure_ROH_compact_total_FROH_summary.tsv"
    )
)

# ============================================================
# PLOT
# ============================================================

p <- ggplot() +

    # --------------------------------------------------------
    # Mean ROH composition
    # --------------------------------------------------------

    geom_col(
        data = class_mean,
        aes(
            x = state,
            y = mean_FROH,
            fill = roh_class
        ),
        width = 0.62,
        colour = NA
    ) +

    # --------------------------------------------------------
    # Seven simulation-replicate values
    # --------------------------------------------------------

    geom_point(
        data = rep_total,
        aes(
            x = state,
            y = total_FROH,
            colour = state
        ),
        position = position_jitter(
            width = 0.09,
            height = 0,
            seed = 1
        ),
        size = 1.8,
        alpha = 0.65
    ) +

    # --------------------------------------------------------
    # SD across simulation replicates
    # --------------------------------------------------------

    geom_errorbar(
        data = total_summary,
        aes(
            x = state,
            ymin = pmax(
                0,
                mean_total - sd_total
            ),
            ymax = mean_total + sd_total,
            colour = state
        ),
        width = 0.13,
        linewidth = 0.65
    ) +

    # --------------------------------------------------------
    # Mean total FROH
    # --------------------------------------------------------

    geom_point(
        data = total_summary,
        aes(
            x = state,
            y = mean_total,
            colour = state
        ),
        shape = 21,
        fill = "white",
        size = 3.0,
        stroke = 1.0
    ) +

    # --------------------------------------------------------
    # Regional panels
    # --------------------------------------------------------

    facet_wrap(
        ~ region,
        nrow = 1
    ) +

    # --------------------------------------------------------
    # Scales
    # --------------------------------------------------------

    scale_fill_manual(
        values = roh_colours,
        drop = FALSE
    ) +

    scale_colour_manual(
        values = scenario_colours,
        guide = "none"
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0, 0.06)
        )
    ) +

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    labs(
        x = NULL,
        y = expression(F[ROH]),
        fill = "ROH length"
    ) +

    # ========================================================
    # NOLEN-STYLE CLEAN THEME
    # ========================================================

    theme_classic(
        base_size = 11,
        base_family = "sans"
    ) +

    theme(

        # No background grid
        panel.grid = element_blank(),

        # Full rectangular box around every population
        panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = 0.55
        ),

        # Region labels in simple boxed strips
        strip.background = element_rect(
            fill = "white",
            colour = "black",
            linewidth = 0.55
        ),

        strip.text = element_text(
            family = "sans",
            face = "bold",
            size = 11,
            margin = margin(
                t = 5,
                r = 4,
                b = 5,
                l = 4
            )
        ),

        # Axis text
        axis.text = element_text(
            family = "sans",
            colour = "black"
        ),

        axis.text.x = element_text(
            size = 8.3,
            margin = margin(
                t = 4
            )
        ),

        axis.text.y = element_text(
            size = 9
        ),

        axis.title.y = element_text(
            family = "sans",
            size = 11,
            margin = margin(
                r = 7
            )
        ),

        # Ticks
        axis.ticks = element_line(
            colour = "black",
            linewidth = 0.4
        ),

        axis.ticks.length = unit(
            2.5,
            "pt"
        ),

        # Legend
        legend.position = "bottom",

        legend.title = element_text(
            family = "sans",
            size = 10
        ),

        legend.text = element_text(
            family = "sans",
            size = 9.5
        ),

        legend.key.height = unit(
            0.45,
            "cm"
        ),

        legend.key.width = unit(
            0.55,
            "cm"
        ),

        # Spacing between regional boxes
        panel.spacing = unit(
            0.65,
            "lines"
        ),

        # Clean margins
        plot.margin = margin(
            t = 6,
            r = 8,
            b = 4,
            l = 5
        )
    )

# ============================================================
# SAVE
# ============================================================

ggsave(
    file.path(
        OUT,
        "Figure_ROH_compact_by_region.svg"
    ),
    p,
    width = 11.3,
    height = 4.0,
    units = "in",
    bg = "white"
)

ggsave(
    file.path(
        OUT,
        "Figure_ROH_compact_by_region.png"
    ),
    p,
    width = 11.3,
    height = 4.0,
    units = "in",
    dpi = 400,
    bg = "white"
)

cat("\n============================================\n")
cat("UPDATED COMPACT ROH FIGURE COMPLETE\n")
cat("============================================\n\n")

cat(
    file.path(
        OUT,
        "Figure_ROH_compact_by_region.svg"
    ),
    "\n"
)
