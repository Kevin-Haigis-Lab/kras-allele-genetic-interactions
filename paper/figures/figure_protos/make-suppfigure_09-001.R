# Supplemental Figure 9. Labeled high-level view of comutation network for MM

FIGNUM <- 9
SUPPLEMENTAL <- TRUE
VERSION <- 1

FIG_DIMENSIONS <- get_figure_dimensions(2, "medium")


#' General theme for Supp. Figure 9.
theme_figS9 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm"))
    )
}


#' Special theme for graphs from 'ggraph'.
theme_graph_figS9 <- function() {
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

#### ---- A. Labeled MM comutation graph ---- ####
# The comutation graph for MM with every node labeled.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_labeled_MM",
                          FIGNUM,
                          supp = SUPPLEMENTAL) +
    scale_edge_color_manual(
        values = comut_updown_pal,
        guide = guide_legend(
            title = "comutation",
            title.position = "top",
            keywidth = unit(2, "mm"),
            keyheight = unit(1, "mm"),
            nrow = 1,
            label.position = "top",
            order = 1,
            direction = "horizontal"
        )
    ) +
    scale_color_manual(
        values = short_allele_pal,
        na.value = NA,
        guide = guide_legend(
            title = NULL,
            keywidth = unit(2, "mm"),
            keyheight = unit(3, "mm"),
            nrow = 2,
            order = 2
        )
    ) +
    theme_graph_figS9() +
    theme(
        legend.position = c(0, 0.1),
        legend.direction = "horizontal",
        legend.justification = "left"
    ) +
    labs(tag = "a")


#### ---- B. Labeled MM comutation graph ---- ####
# The comutation graph for MM with every node labeled.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_B_1 <- read_fig_proto("mm_comut_heatmap", FIGNUM, supp = SUPPLEMENTAL) +
    theme_figS9() +
    theme(
        legend.position = "bottom",
        axis.title = element_blank()
    )
panel_B_2 <- read_fig_proto("allele_freq_barplot", FIGNUM, supp = SUPPLEMENTAL) +
    theme_figS9() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
    ) +
    labs(tag = "b")
panel_B_3 <- read_fig_proto("gene_freq_barplot", FIGNUM, supp = SUPPLEMENTAL) +
    theme_figS9() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
    )

panel_B <-
    (panel_B_2 + plot_spacer() + plot_layout(widths = c(10, 2))) /
    (panel_B_1 + panel_B_3 + plot_layout(widths = c(10, 2))) +
    plot_layout(heights = c(2, 10))

#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- panel_A / (panel_B) +
        plot_layout(heights = c(2, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
