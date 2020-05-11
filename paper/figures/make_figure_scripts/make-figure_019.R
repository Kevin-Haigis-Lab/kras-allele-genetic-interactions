# Figure 019. Analysis of the JNK and c-Jun dependencies in PAAD.

FIGNUM <- 19

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")
FIG_DIMENSIONS$height <- FIG_DIMENSIONS$height / 2


theme_fig19 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. Gene set enrichment heatmap for JNK activation ---- ####
# A heatmap of an enriched gene set with the cell lines ranked by their
# dependency score for each gene. A density plot along the top helps
# highlight the trend.
# original script: "src/10_37_gsea-depmap-analysis.R"

theme_fig19_densityplots <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_classic_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin
        ),
        plot.title = element_text(size = 6, family = "Arial", vjust = 1),
        plot.margin = margin(0, 0, -1.5, 0, "mm"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
    )
}

# The proto object is actually in the directory for Supp. Fig 13 with
# the rest of the PAAD files.
PAAD_MAIN_FIGNUM <- 13

x_label <- expression("" %<-% "greater dep. - ranked gene effect - less dep." %->% "")


panel_A_density <- read_fig_proto(
        "rankline_PAAD_G12V_REACTOME_JNK_C_JUN_KINASES_PHOSPHORYLATION_AND_ACTIVATION_MEDIATED_BY_ACTIVATED_HUMAN_TAK1.rds"
    ) +
    theme_fig19_densityplots() +
    labs(
        y = "density",
        title = "JNK phosphorylation and activation by activated TAK1"
    )

panel_A <- read_fig_proto(
        "rankplot_PAAD_G12V_REACTOME_JNK_C_JUN_KINASES_PHOSPHORYLATION_AND_ACTIVATION_MEDIATED_BY_ACTIVATED_HUMAN_TAK1.rds"
    ) +
    theme_fig19() +
    theme(
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = NA, color = NA),
        plot.margin = margin(0, 0, 0, 0, "mm")
    ) +
    labs(
        x = x_label,
        fill = "KRAS allele"
    )


#### ---- B. Genetic dependency box-plot ---- ####
# Box-plot of MAPK8 genetic dependency.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"


panel_B <- read_fig_proto("PAAD-MAPK8.rds") +
    theme_fig19(tag_margin = margin(0, 2, 0, -2, "mm")) %+replace%
    theme(
        plot.title = element_text(size = 6, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 3, 0, "mm")
    ) +
    labs(
        title = "MAPK8",
        tag = "b"
    )


#### ---- C. Genetic dependency scatter plot for JUN vs. MAPK8 ---- ####
# Scatter lpot of JUN and MAPK* genetic dependency.
# original script: "src/10_55_paad_depmap_jun-cdkn2a-G12V.R"


panel_C <- read_fig_proto("JUN_MAPK8_scatter.rds") +
    theme_fig19(tag_margin = margin(0, 2, 0, -2, "mm")) %+replace%
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 3, 0, "mm")
    ) +
    labs(
        x = "dep. on JUN",
        y = "dep. on MAPK8",
        tag = "c"
    )


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE

    col1 <- wrap_elements(
        full = (
                plot_spacer() / panel_A_density / panel_A / plot_spacer()
            ) +
            plot_layout(heights = c(3, 5, 11, 1))
    ) +
        theme_fig19() +
        labs(tag = "a")

    col2 <- wrap_elements(
        full = (panel_B / panel_C)
    )

    full_figure <- (plot_spacer() | col1 | col2 | plot_spacer()) +
        plot_layout(widths = c(1, 9, 4, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
