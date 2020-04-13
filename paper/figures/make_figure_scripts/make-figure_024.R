# Figure 024. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 24

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 90

theme_fig24 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


# Special theme for graphs from 'ggraph'.
theme_graph_fig24 <- function() {
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


#### ---- A. High-level comutation network for MM (labeled) ---- ####
# The labeled, high-level network plot for the comutation graph for MM
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_labeled_MM") +
    theme_graph_fig24() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = c(0.1, 0.1)
    ) +
    labs(title = "MM")



#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (plot_spacer() | panel_A | plot_spacer()) +
        plot_layout(widths = c(1, 3, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
