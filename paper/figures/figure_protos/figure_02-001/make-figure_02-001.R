# Build Figure 2

FIGNUM <- 2
VERSION <- 1
FIGFILENAME <- glue("figure_{FIGNUM}_{VERSION}.svg")
FIG_DIMENSIONS <- get_figure_dimensions(2, "medium")


library(patchwork)


#### ---- Figure theme ---- ####

theme_fig2 <- function() {
    theme_bw(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5, hjust = 1),
        axis.text.x = element_text(size = 5, vjust = 1),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-3, -3, -3, -3, "mm"))
    )
}

theme_graph_fig2 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        legend.position = "left",
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-3, -3, -3, -3, "mm"))
    )
}


#### ---- A. High-level comutation network for COAD ---- ####

# Panel A.
# The high-level network plot for the comutation graph for COAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_COAD.svg", 2) *
    theme_graph_fig2()


#### ---- B. A priori genes of interest comutation network for COAD ---- ####

# Panel B.
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"

panel_B <- read_fig_proto("goi_overlap_genetic_interactions_network_COAD_allLists", 2) *
    theme_graph_fig2() +
    theme(
        legend.position = "bottom"
    )



#### ---- Figure assembly ---- ####

{
    # ROW 1
    row_1 <- panel_A + panel_B + plot_spacer() / patchwork::plot_spacer() +
        plot_layout(heights = c(1, 2))

    # COMPLETE FIGURE
    full_figure <- row_1 +
        plot_annotation(
            title = glue("Figure {FIGNUM}"),
            theme = theme(
                plot.title = element_text(size = 10, family = "Arial")
            )
        )

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
