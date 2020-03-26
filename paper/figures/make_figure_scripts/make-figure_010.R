# Figure 010. Labeled high-level view of comutation network for COAD

FIGNUM <- 10

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


#' Special theme for graphs from 'ggraph'.
theme_graph_fig10 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm")),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "bottom",
        legend.title = element_text(size = 5, hjust = 0.5),
        legend.text = element_text(size = 5, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}


#### ---- A. Labeled COAD comutation graph ---- ####
# The comutation graph for COAD with every node labeled.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_labeled_COAD") +
    theme_graph_fig10()


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- panel_A +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
