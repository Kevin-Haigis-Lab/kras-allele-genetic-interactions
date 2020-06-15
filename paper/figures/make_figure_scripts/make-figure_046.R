# Figure 046. #> Second half of results from explaining allele-specific
# genetic dependency with comutation interactions.

FIGNUM <- 46

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 120

theme_fig46 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}



theme_coefplot_fig46 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_minimal_comutation() %+replace%
    theme(
        plot.title = element_markdown(size = 7,
                                      face = "bold",
                                      family = "Arial",
                                      hjust = 0),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown()
    )
}


theme_lineplot_fig46 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_fig45(tag_margin = tag_margin) %+replace%
    theme(
        legend.position = "bottom",
        axis.title.y = element_markdown(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.key.width = unit(2, "mm")
    )
}



theme_fig46_void <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_graph_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


subpanel_heights <- c(2, 8, 2)


break_labels_ampersand <- function(x) {
    str_replace(x, " \\& ", " &\n")
}


#### ---- A. COAD - G12D - SRSF5 ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("COAD_G12D_SRSF5_coef-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.15))) +
    scale_y_discrete(expand = expansion(mult = c(0.3, 0.3))) +
    theme_coefplot_fig46() +
    labs(tag = "a")
panel_A_2 <- read_fig_proto("COAD_G12D_SRSF5_box-plot") +
    scale_x_discrete(labels = break_labels_ampersand) +
    theme_lineplot_fig46()
panel_A_3 <- read_fig_proto("COAD_G12D_SRSF5_comut-graph-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    theme_fig46_void()


panel_A <- (panel_A_1 / panel_A_2 / panel_A_3) +
    plot_layout(heights = subpanel_heights)



#### ---- B. PAAD - G12R - KIAA1257 ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12R_KIAA1257_coef-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.05))) +
    theme_coefplot_fig46() +
    labs(tag = "b")
panel_B_2 <- read_fig_proto("PAAD_G12R_KIAA1257_box-plot") +
    scale_x_discrete(labels = break_labels_ampersand) +
    theme_lineplot_fig46() +
    theme(
    )
panel_B_3 <- read_fig_proto("PAAD_G12R_KIAA1257_comut-graph-plot.rds") +
    theme_fig46_void()

panel_B <- panel_B_1 / panel_B_2 / panel_B_3 +
    plot_layout(heights = subpanel_heights)



#### ---- C. PAAD - G12D - FKBP1A ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_C_1 <- read_fig_proto("PAAD_G12D_FKBP1A_coef-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.25, 0.05))) +
    theme_coefplot_fig46() +
    labs(tag = "c")
panel_C_2 <- read_fig_proto("PAAD_G12D_FKBP1A_box-plot") +
    scale_x_discrete(labels = break_labels_ampersand, drop = FALSE) +
    theme_lineplot_fig46()
panel_C_3 <- read_fig_proto("PAAD_G12D_FKBP1A_comut-graph-plot.rds") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    theme_fig46_void()

panel_C <- panel_C_1 / panel_C_2 / panel_C_3 +
    plot_layout(heights = subpanel_heights)



#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (panel_A | panel_B | panel_C)

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
