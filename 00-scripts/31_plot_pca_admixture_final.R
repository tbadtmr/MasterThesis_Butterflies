#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(patchwork)
})

# ============================================================
# PATHS
# ============================================================

BASE <- path.expand(
    "~/tabea_work/08_population_analysis_final/skane_only/combined/results/structure/final_plot"
)

PCA_FILE <- file.path(
    BASE,
    "representative_PCA.tsv"
)

Q_FILE <- file.path(
    BASE,
    "representative_admixture_Q.tsv"
)

K_FILE <- file.path(
    BASE,
    "admixture_K_summary.tsv"
)

SEL_FILE <- file.path(
    BASE,
    "admixture_selection.txt"
)

OUT <- BASE

pca <- read_tsv(
    PCA_FILE,
    show_col_types = FALSE
)

q <- read_tsv(
    Q_FILE,
    show_col_types = FALSE
)

ks <- read_tsv(
    K_FILE,
    show_col_types = FALSE
)

# ============================================================
# ORDERS
# ============================================================

state_order <- c(
    "1900",
    "2020",
    "SQ",
    "R2",
    "R4",
    "R6"
)

state_titles <- c(
    "1900",
    "2020",
    "2140 Status quo",
    "2140 Restore 2 km",
    "2140 Restore 4 km",
    "2140 Restore 6 km"
)

names(state_titles) <- state_order

region_order <- c(
    "NW",
    "NE",
    "W",
    "E",
    "SE"
)

# ============================================================
# YOUR SCENARIO COLOURS
# ============================================================

scenario_colours <- c(
    "1900" = "#8A8A8A",
    "2020" = "#555555",
    "SQ"   = "#000000",
    "R2"   = "#69A541",
    "R4"   = "#21A6B6",
    "R6"   = "#9A6DCE"
)

# ============================================================
# REGIONAL COLOURS FOR PCA
#
# Restrained geographic palette, inspired by the empirical
# structure figure rather than by scenario colours.
# ============================================================

region_colours <- c(
    "NW" = "#244A6A",
    "NE" = "#527C9B",
    "W"  = "#7FA5B8",
    "E"  = "#A8616A",
    "SE" = "#7D3044"
)

# ============================================================
# ADMIXTURE CLUSTER COLOURS
#
# Kept separate from scenario colours.
# Supports K up to 8.
# ============================================================

cluster_colours <- c(
    "#244A6A",
    "#8FB3C9",
    "#84384A",
    "#D5A72D",
    "#737373",
    "#6F8F63",
    "#8467A9",
    "#3F8B8B"
)

# ============================================================
# COMMON DATA PREP
# ============================================================

pca <- pca %>%
    mutate(
        state = factor(
            state,
            levels = state_order
        ),

        region = factor(
            region,
            levels = region_order
        )
    )

q <- q %>%
    mutate(
        state = factor(
            state,
            levels = state_order
        ),

        region = factor(
            region,
            levels = region_order
        )
    )

q_cols <- grep(
    "^Q[0-9]+$",
    names(q),
    value = TRUE
)

BEST_K <- length(q_cols)

cat(
    "Plotting ADMIXTURE K =",
    BEST_K,
    "\n"
)

# ============================================================
# COMMON PCA AXIS LIMITS
# ============================================================

pc1_range <- range(pca$PC1, na.rm = TRUE)
pc2_range <- range(pca$PC2, na.rm = TRUE)

pc1_pad <- diff(pc1_range) * 0.04
pc2_pad <- diff(pc2_range) * 0.04

pca_xlim <- c(
    pc1_range[1] - pc1_pad,
    pc1_range[2] + pc1_pad
)

pca_ylim <- c(
    pc2_range[1] - pc2_pad,
    pc2_range[2] + pc2_pad
)

cat(
    "Common PCA PC1 limits:",
    paste(round(pca_xlim, 4), collapse = " to "),
    "\n"
)

cat(
    "Common PCA PC2 limits:",
    paste(round(pca_ylim, 4), collapse = " to "),
    "\n"
)

# ============================================================
# A. PCA
# ============================================================

make_pca <- function(st) {

    d <- pca %>%
        filter(
            state == st
        )

    ggplot(
        d,
        aes(
            x = PC1,
            y = PC2,
            colour = region
        )
    ) +

        geom_point(
            size = 1.5,
            alpha = 0.72
        ) +

        scale_colour_manual(
            values = region_colours,
            drop = FALSE
        ) +

        scale_x_continuous(
            limits = pca_xlim,
            breaks = scales::breaks_pretty(n = 4)
        ) +

        scale_y_continuous(
            limits = pca_ylim,
            breaks = scales::breaks_pretty(n = 4)
        ) +

        labs(
            title = state_titles[[st]],
            x = "PC1",
            y = "PC2",
            colour = "Region",
            tag = if (st == "1900") "A" else NULL
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
                colour = scenario_colours[[st]],
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

            legend.position = "none",

            plot.margin = margin(
                3, 3, 3, 3
            )
        )
}

pca_plots <- lapply(
    state_order,
    make_pca
)

pca_row <- wrap_plots(
    pca_plots,
    nrow = 1
)

# ============================================================
# B. ADMIXTURE PREPARATION
#
# Individuals ordered by region, site, then individual index.
# ============================================================

q <- q %>%
    arrange(
        state,
        region,
        site_number,
        site,
        ind_index
    ) %>%

    group_by(
        state
    ) %>%

    mutate(
        plot_x = row_number()
    ) %>%

    ungroup()

qlong <- q %>%

    pivot_longer(
        cols = all_of(q_cols),
        names_to = "ancestry_cluster",
        values_to = "ancestry"
    ) %>%

    mutate(
        ancestry_cluster = factor(
            ancestry_cluster,
            levels = q_cols
        )
    )

# ============================================================
# REGION BOUNDARIES AND LABEL POSITIONS
# ============================================================

region_positions <- q %>%

    group_by(
        state,
        region
    ) %>%

    summarise(
        xmin = min(plot_x),
        xmax = max(plot_x),
        midpoint = (
            min(plot_x) +
            max(plot_x)
        ) / 2,
        .groups = "drop"
    )

# ============================================================
# ADMIXTURE PLOT FUNCTION
# ============================================================

make_admix <- function(st) {

    d <- qlong %>%
        filter(
            state == st
        )

    pos <- region_positions %>%
        filter(
            state == st
        )

    ggplot(
        d,
        aes(
            x = plot_x,
            y = ancestry,
            fill = ancestry_cluster
        )
    ) +

        geom_col(
            width = 1,
            colour = NA
        ) +

        # vertical boundaries between regional populations
        geom_vline(
            data = pos[-nrow(pos), ],
            aes(
                xintercept = xmax + 0.5
            ),
            linewidth = 0.35,
            colour = "white"
        ) +

        # regional labels above bars
        geom_text(
            data = pos,
            aes(
                x = midpoint,
                y = 1.055,
                label = region,
                colour = region
            ),
            inherit.aes = FALSE,
            size = 2.7,
            fontface = "bold",
            family = "sans"
        ) +

        scale_colour_manual(
            values = region_colours,
            drop = FALSE
        ) +

        scale_fill_manual(
            values = cluster_colours[
                seq_len(BEST_K)
            ],
            drop = FALSE
        ) +

        scale_y_continuous(
            limits = c(
                0,
                1.09
            ),
            breaks = c(
                0,
                0.5,
                1
            ),
            expand = c(
                0,
                0
            )
        ) +

        scale_x_continuous(
            expand = c(
                0,
                0
            )
        ) +

        labs(
            title = state_titles[[st]],
            x = NULL,
            y = "Ancestry",
            fill = NULL,
            colour = "Region",
            tag = if (st == "1900") "B" else NULL
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
                colour = scenario_colours[[st]],
                face = "bold",
                size = 9.5,
                hjust = 0.5,
                margin = margin(
                    b = 4
                )
            ),

            axis.text.x =
                element_blank(),

            axis.ticks.x =
                element_blank(),

            axis.title.y =
                element_text(
                    size = 8
                ),

            axis.text.y =
                element_text(
                    size = 7
                ),

            legend.position =
                "none",

            plot.margin =
                margin(
                    3, 3, 3, 3
                )
        )
}

admix_plots <- lapply(
    state_order,
    make_admix
)

admix_row <- wrap_plots(
    admix_plots,
    nrow = 1
)

# ============================================================
# CLUSTER LEGEND
# ============================================================

legend_dat <- data.frame(
    cluster = factor(
        q_cols,
        levels = q_cols
    ),
    x = seq_along(q_cols),
    y = 1
)

cluster_legend <- ggplot(
    legend_dat,
    aes(
        x = x,
        y = y,
        fill = cluster
    )
) +

    geom_tile() +

    scale_fill_manual(
        values = cluster_colours[
            seq_len(BEST_K)
        ],
        labels = paste0(
            "Cluster ",
            seq_len(BEST_K)
        )
    ) +

    guides(
        fill = guide_legend(
            nrow = 1,
            title = NULL
        )
    ) +

    theme_void() +

    theme(
        legend.position = "bottom",
        legend.text = element_text(
            family = "sans",
            size = 8
        )
    )

# ============================================================
# FINAL FIGURE
# ============================================================

final_plot <- (
    pca_row /
    admix_row
) +

    plot_layout(
        heights = c(
            1.35,
            0.80
        ),
        guides = "collect"
    ) +

    plot_annotation() &

    theme(
        legend.position = "bottom"
    )

# ============================================================
# SAVE
# ============================================================

ggsave(
    file.path(
        OUT,
        "Figure_PCA_ADMIXTURE_final.svg"
    ),
    final_plot,
    width = 12.5,
    height = 6.2,
    units = "in",
    bg = "white"
)

ggsave(
    file.path(
        OUT,
        "Figure_PCA_ADMIXTURE_final.png"
    ),
    final_plot,
    width = 12.5,
    height = 6.2,
    units = "in",
    dpi = 400,
    bg = "white"
)

# ============================================================
# CV K-SELECTION FIGURE
# ============================================================

best_k <- ks$K[
    which.min(
        ks$mean_cv
    )
]

p_cv <- ggplot(
    ks,
    aes(
        x = K,
        y = mean_cv
    )
) +

    geom_errorbar(
        aes(
            ymin = mean_cv - sd_cv,
            ymax = mean_cv + sd_cv
        ),
        width = 0.15,
        linewidth = 0.45
    ) +

    geom_line(
        linewidth = 0.55
    ) +

    geom_point(
        size = 2
    ) +

    geom_point(
        data = ks %>%
            filter(
                K == best_k
            ),
        size = 3.5,
        shape = 21,
        fill = "white",
        stroke = 1
    ) +

    scale_x_continuous(
        breaks = ks$K
    ) +

    labs(
        x = "Number of ancestral clusters (K)",
        y = "Cross-validation error"
    ) +

    theme_classic(
        base_family = "sans",
        base_size = 10
    ) +

    theme(
        panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = 0.45
        )
    )

ggsave(
    file.path(
        OUT,
        "Figure_ADMIXTURE_K_selection.svg"
    ),
    p_cv,
    width = 5,
    height = 3.6,
    units = "in",
    bg = "white"
)

ggsave(
    file.path(
        OUT,
        "Figure_ADMIXTURE_K_selection.png"
    ),
    p_cv,
    width = 5,
    height = 3.6,
    units = "in",
    dpi = 400,
    bg = "white"
)

cat("\n============================================\n")
cat("PCA + ADMIXTURE FIGURES COMPLETE\n")
cat("============================================\n\n")

cat(
    "Best K:",
    best_k,
    "\n"
)

cat(
    "Final figure:\n",
    file.path(
        OUT,
        "Figure_PCA_ADMIXTURE_final.svg"
    ),
    "\n"
)

cat(
    "K-selection figure:\n",
    file.path(
        OUT,
        "Figure_ADMIXTURE_K_selection.svg"
    ),
    "\n"
)
