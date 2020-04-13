# Figure 023. #> LUAD comutation network

FIGNUM <- 23

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 90

theme_fig23 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


# Special theme for graphs from 'ggraph'.
theme_graph_fig23 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, face = "bold",
                                  hjust = 0.1, vjust = 5),
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


#### ---- A. High-level comutation network for LUAD ---- ####
# The high-level network plot for the comutation graph for LUAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_LUAD") +
    theme_graph_fig23() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = c(0.05, 0.00),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6)
    ) +
    labs(tag = "a",
         title = "LUAD")



#### ---- B. A priori genes of interest comutation network for LUAD ---- ####
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"
panel_B <- read_fig_proto(
        "goi_overlap_genetic_interactions_network_LUAD_kegg"
    ) +
    scale_size_manual(
        values = c(big = 1.6, small = 1.5),
        guide = FALSE
    ) +
    theme_graph_fig23() %+replace%
    theme(
        legend.title = element_markdown(size = 6),
        legend.text = element_text(size = 6)
    ) +
    labs(
        tag = "b",
        edge_width = "-*log*(p-value)",
        title = "LUAD"
    )


#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    # COMPLETE FIGURE
    full_figure <- (panel_A | panel_B)

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
