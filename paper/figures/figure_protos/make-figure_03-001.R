# Build Figure 3

FIGNUM <- 3
VERSION <- 1
FIGFILENAME <- glue("figure_{FIGNUM}_{VERSION}.svg")
FIG_DIMENSIONS <- get_figure_dimensions(1, "tall")

library(ggraph)


#### ---- Figure theme ---- ####

#' General theme for Figure 2.
theme_fig3 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm"))
    )
}


#' Special theme for graphs from 'ggraph'.
theme_graph_fig3 <- function() {
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



#### ---- A. High-level comutation network for LUAD ---- ####

# Panel A.
# The high-level network plot for the comutation graph for LUAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_LUAD", FIGNUM) +
    theme_graph_fig3() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm")
    ) +
    labs(tag = "a")


#### ---- B. Dot-plot of functional enrichment ---- ####

# Panel B.
# A dot plot of the results of functional enrichment in LUAD comutation network.
# original script: "src/20_45_fxnal-enrich-genetic-interactions.R"
panel_B <- read_fig_proto("enrichr_LUAD", FIGNUM) +
    scale_size_continuous(
        range = c(0, 3),
        guide = guide_legend(
            title.position = "left",
            title.hjust = 0.5,
            order = 10,
            label.vjust = -10,
            label.position = "top"
        )
    ) +
    scale_alpha_continuous(
        range = c(0, 1),
        guide = guide_legend(
            title.position = "left",
            title.hjust = 0.5,
            order = 20,
            label.vjust = -10,
            label.position = "top"
        )
    ) +
    theme_fig2() +
    theme(
        plot.title = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(-1, 0, -1, 0, "mm"),
        legend.box.background = element_rect(fill = NA, color = NA)
    ) +
    labs(tag = "c")


#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    # COMPLETE FIGURE
    full_figure <- (
        panel_A /
        ((plot_spacer() + panel_B) + plot_layout(widths = c(1, 1e6))) /
        (plot_spacer() + plot_spacer() + plot_layout(widths = c(1, 1e6)))
    ) +
        plot_layout(heights = c(1, 1, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}


