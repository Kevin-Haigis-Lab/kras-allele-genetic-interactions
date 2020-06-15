# Figure 045. #> The first group of results from the analysis for explaining
# allele-specific genetic dependency by comutation.

FIGNUM <- 45

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 120


theme_fig45 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


theme_coefplot_fig45 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
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


theme_lineplot_fig45 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
    theme_fig45(tag_margin = tag_margin) %+replace%
    theme(
        legend.position = "none",
        axis.title.y = element_markdown(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_blank()
    )
}


theme_fig45_void <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_graph_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. TP53 explaining STARD9 ---- ####
# TP53 explaining dep. of G12D on STARD9 in COAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("COAD_G12D_STARD9_coef-plot") +
    theme_coefplot_fig45() +
    labs(tag = "a")
panel_A_2 <- read_fig_proto("COAD_G12D_STARD9_line-plot") +
    theme_lineplot_fig45()
panel_A_3 <- read_fig_proto("COAD_G12D_TP53_comut-graph-plot.rds") +
    scale_x_continuous(expand = expansion(mult = c(0.50, 0.50))) +
    theme_fig45_void()



#### ---- B. SMAD4 explaining EEF1E1 ---- ####
# SMAD4 explaining dep. of G12D on EEF1E1 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12D_EEF1E1_coef-plot") +
    theme_coefplot_fig45() +
    labs(tag = "b")
panel_B_2 <- read_fig_proto("PAAD_G12D_EEF1E1_line-plot") +
    theme_lineplot_fig45()
panel_B_3 <- read_fig_proto("PAAD_G12D_SMAD4_comut-graph-plot.rds") +
    scale_x_continuous(expand = expansion(mult = c(0.50, 0.50))) +
    theme_fig45_void()


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (
        (panel_A_1 / panel_A_2 / panel_A_3) + plot_layout(height = c(2, 8, 2)) |
        (panel_B_1 / panel_B_2 / panel_B_3) + plot_layout(height = c(2, 8, 2))
    )

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
