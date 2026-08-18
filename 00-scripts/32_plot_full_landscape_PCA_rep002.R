#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# ============================================================
# PATHS
# ============================================================

BASE <- path.expand(
    "~/tabea_work/08_population_analysis"
)

INFILE <- file.path(
    BASE,
    "plots",
    "k13_full",
    "regional_pca_joint",
    "regional20",
    "ALL_full_joint_PCA_scores.tsv"
)

OUTDIR <- file.path(
    BASE,
    "plots",
    "k13_full",
    "regional_pca_joint",
    "regional20",
    "supplement_final"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(INFILE)) {
    stop(
        "Input file not found:\n",
        INFILE
    )
}

# ============================================================
# SETTINGS
# ============================================================

REP_TO_PLOT <- "rep002"

state_order <- c(
    "1900",
    "2020",
    "2140 Status quo",
    "2140 2 km restoration",
    "2140 4 km restoration",
    "2140 6 km restoration"
)

state_titles <- c(
    "1900"                  = "1900",
    "2020"                  = "2020",
    "2140 Status quo"       = "2140 Status quo",
    "2140 2 km restoration" = "2140 Restore 2 km",
    "2140 4 km restoration" = "2140 Restore 4 km",
    "2140 6 km restoration" = "2140 Restore 6 km"
)

# Same scenario-title colours as the new Skåne figure
scenario_colours <- c(
    "1900"                  = "#8A8A8A",
    "2020"                  = "#555555",
    "2140 Status quo"       = "#000000",
    "2140 2 km restoration" = "#69A541",
    "2140 4 km restoration" = "#21A6B6",
    "2140 6 km restoration" = "#9A6DCE"
)

# ============================================================
# REGION LABELS
# ============================================================

region_label_map <- c(
    "NW"    = "NW",
    "NE"    = "NE",
    "W"     = "W",
    "E"     = "E",
    "SE"    = "SE",
    "NEW01" = "Öland",
    "NEW02" = "W Småland"
)

region_order <- c(
    "NW",
    "NE",
    "W",
    "E",
    "SE",
    "W Småland",
    "Öland"
)

# Same five Skåne colours as the new PCA + two additional
# colours for W Småland and Öland.
region_colours <- c(
    "NW"        = "#244A6A",
    "NE"        = "#527C9B",
    "W"         = "#7FA5B8",
    "E"         = "#A8616A",
    "SE"        = "#7D3044",
    "W Småland" = "#666666",
    "Öland"     = "#7B3294"
)

# ============================================================
# READ DATA
# ============================================================

dat <- read_tsv(
    INFILE,
    show_col_types = FALSE
) %>%

    filter(
        rep == REP_TO_PLOT
    ) %>%

    mutate(
        region_plot = recode(
            region,
            !!!region_label_map
        ),

        region_plot = factor(
            region_plot,
            levels = region_order
        ),

        state = factor(
            state,
            levels = state_order
        )
    )

if (nrow(dat) == 0) {
    stop(
        "No data found for ",
        REP_TO_PLOT
    )
}

cat(
    "Individuals plotted:",
    nrow(dat),
    "\n"
)

cat(
    "Regions:",
    paste(
        unique(as.character(dat$region_plot)),
        collapse = ", "
    ),
    "\n"
)

# ============================================================
# COMMON PCA LIMITS
#
# Important: all six states came from ONE joint PCA.
# Therefore they must use the same plotting axes.
# ============================================================

pc1_range <- range(
    dat$PC1,
    na.rm = TRUE
)

pc2_range <- range(
    dat$PC2,
    na.rm = TRUE
)

pc1_pad <- diff(pc1_range) * 0.04
pc2_pad <- diff(pc2_range) * 0.04

if (pc1_pad == 0) pc1_pad <- 0.01
if (pc2_pad == 0) pc2_pad <- 0.01

pca_xlim <- c(
    pc1_range[1] - pc1_pad,
    pc1_range[2] + pc1_pad
)

pca_ylim <- c(
    pc2_range[1] - pc2_pad,
    pc2_range[2] + pc2_pad
)

cat(
    "Common PC1 limits:",
    paste(
        round(pca_xlim, 4),
        collapse = " to "
    ),
    "\n"
)

cat(
    "Common PC2 limits:",
    paste(
        round(pca_ylim, 4),
        collapse = " to "
    ),
    "\n"
)

# ============================================================
# AXIS LABELS
#
# These are the variance percentages from the existing
# rep002 six-state joint PCA figure.
# ============================================================

PC1_LABEL <- "PC1 (68.9%)"
PC2_LABEL <- "PC2 (9.0%)"

# ============================================================
# PLOT ONE STATE
# ============================================================

make_pca <- function(st) {

    d <- dat %>%
        filter(
            state == st
        )

    ggplot(
        d,
        aes(
            x = PC1,
            y = PC2,
            colour = region_plot
        )
    ) +

        geom_point(
            size = 1.55,
            alpha = 0.72
        ) +

        scale_colour_manual(
            values = region_colours,
            breaks = region_order,
            drop = FALSE,
            name = "Region"
        ) +

        scale_x_continuous(
            limits = pca_xlim,
            breaks = scales::breaks_pretty(
                n = 4
            )
        ) +

        scale_y_continuous(
            limits = pca_ylim,
            breaks = scales::breaks_pretty(
                n = 4
            )
        ) +

        labs(
            title = state_titles[[st]],
            x = PC1_LABEL,
            y = PC2_LABEL
        ) +

        theme_classic(
            base_size = 9,
            base_family = "sans"
        ) +

        theme(
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 0.45
            ),

            plot.title = element_text(
                colour =
                    scenario_colours[[st]],
                face = "bold",
                size = 9.5,
                hjust = 0.5,
                margin = margin(
                    b = 4
                )
            ),

            axis.title = element_text(
                size = 8
            ),

            axis.text = element_text(
                size = 7,
                colour = "black"
            ),

            legend.position = "bottom",

            plot.margin = margin(
                3,
                3,
                3,
                3
            )
        )
}

# ============================================================
# BUILD SIX PANELS
# ============================================================

plots <- lapply(
    state_order,
    make_pca
)

p_row <- wrap_plots(
    plots,
    nrow = 1
)

# ============================================================
# CREATE ONE CLEAN SHARED LEGEND
# ============================================================

legend_dat <- data.frame(
    PC1 = seq_along(region_order),
    PC2 = 1,
    region_plot = factor(
        region_order,
        levels = region_order
    )
)

legend_plot <- ggplot(
    legend_dat,
    aes(
        PC1,
        PC2,
        colour = region_plot
    )
) +

    geom_point(
        size = 2
    ) +

    scale_colour_manual(
        values = region_colours,
        breaks = region_order,
        drop = FALSE,
        name = "Region"
    ) +

    guides(
        colour = guide_legend(
            nrow = 1,
            byrow = TRUE,
            override.aes = list(
                size = 2.5,
                alpha = 1
            )
        )
    ) +

    theme_void(
        base_family = "sans"
    ) +

    theme(
        legend.position = "bottom",
        legend.title = element_text(
            size = 9
        ),
        legend.text = element_text(
            size = 8
        )
    )

legend <- cowplot::get_legend(
    legend_plot
)

# ============================================================
# FINAL FIGURE
#
# All panels use the same region scale, so patchwork collects
# them into one shared legend.
# ============================================================

final_plot <- wrap_plots(
    plots,
    nrow = 1,
    guides = "collect"
) &

    theme(
        legend.position = "bottom",
        legend.title = element_text(
            size = 9
        ),
        legend.text = element_text(
            size = 8
        )
    )

# ============================================================
# SAVE
# ============================================================

SVG <- file.path(
    OUTDIR,
    "Figure_S_full_landscape_PCA_rep002.svg"
)

PNG <- file.path(
    OUTDIR,
    "Figure_S_full_landscape_PCA_rep002.png"
)

ggsave(
    SVG,
    final_plot,
    width = 12.5,
    height = 3.65,
    units = "in",
    bg = "white"
)

ggsave(
    PNG,
    final_plot,
    width = 12.5,
    height = 3.65,
    units = "in",
    dpi = 400,
    bg = "white"
)

cat("\n============================================\n")
cat("FULL-LANDSCAPE PCA COMPLETE\n")
cat("============================================\n\n")

cat(
    "SVG:\n",
    SVG,
    "\n\n"
)

cat(
    "PNG:\n",
    PNG,
    "\n"
)
