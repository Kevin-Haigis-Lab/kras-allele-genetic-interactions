# Figure 031. Heatmap of dependency analysis results for LUAD.

FIGNUM <- 31

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(1, "tall")


theme_fig31 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. Heatmap of linear model ---- ####
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal in LUAD.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_A <- read_fig_proto(
    "LUAD_CRISPR_manhattan_ward.D2_pheatmap.rds"
)[[4]]

pre_panel_A_main <- gtable::gtable_filter(pre_panel_A,
                                          "legend",
                                          invert = TRUE)

panel_A <- wrap_elements(plot = pre_panel_A_main) *
    theme_fig31() +
    theme(
        plot.margin = margin(-12, -9, -15, -8, "mm")
    )

panel_A_legend1 <- prep_pheatmap_colorbar(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap_heatpal.rds"
    ) +
    labs(title = "scaled\ndep. score")

panel_A_legend2 <- prep_pheatmap_legend(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap_allelepal.rds"
    ) +
    labs(title = "allele")

panel_A_legend3 <- prep_pheatmap_legend(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap_clusterpal.rds"
    ) +
    labs(title = "cluster")


panel_A_legend <- (
        panel_A_legend1 / panel_A_legend2 / panel_A_legend3 / plot_spacer()
    ) +
        plot_layout(heights = c(1, 1, 1, 3))

panel_A_legend <- wrap_elements(plot = panel_A_legend) +
    theme_fig31() %+replace%
    theme(
        plot.margin = margin(0, 0, 0, 0, "mm")
    )


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (
            panel_A + (
                (plot_spacer() / panel_A_legend / plot_spacer()) +
                plot_layout(heights = c(1, 2, 1))
            )
        ) + 
        plot_layout(widths = c(5, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
