# Figure 032. PPI subnetwork from LUAD DepMap analysis.

FIGNUM <- 32

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 90

theme_fig32 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. PPI of genes from cluster in heatmap ---- ####
# A large PPI composed of genes in cluster 4 of the heatmap (panel h).
# original script: "src/10_15_linear-modeling-syn-let_ppi-subnetworks.R"

panel_A <- read_fig_proto("LUAD_cluster-4_component-1") +
    theme_graph_comutation()


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- panel_A +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
