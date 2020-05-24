# Figure 046. #> Second half of results from explaining allele-specific
# genetic dependency with comutation interactions.

FIGNUM <- 46

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig46 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. COAD - G12D - SRSF5 ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("COAD_G12D_SRSF5_coef-plot")
panel_A_2 <- read_fig_proto("COAD_G12D_SRSF5_box-plot")
panel_A_3 <- read_fig_proto("COAD_G12D_SRSF5_comut-graph-plot")


panel_A <- (srsf5_coef_plot / srsf5_box_plot / srsf5_comut_plot) +
    plot_layout(heights = c(3, 9, 1))



#### ---- B. PAAD - G12R - KIAA1257 ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12R_KIAA1257_coef-plot")
panel_B_2 <- read_fig_proto("PAAD_G12R_KIAA1257_box-plot")
panel_B_3 <- read_fig_proto("PAAD_G12R_KIAA1257_G12R-DNAH5_line-plot")


panel_B_layout <- c(
    area(t = 1, l = 1, b = 4, r = 10),
    area(t = 5, l = 1, b = 14, r = 10),
    area(t = 6, l = 5.5, b = 7, r = 10)
)

panel_B <- panel_B_1 + panel_B_2 + panel_B_3 +
    plot_layout(design = panel_B_layout)


#### ---- C. PAAD - G12D - FKBP1A ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_C_1 <- read_fig_proto("PAAD_G12D_FKBP1A_coef-plot")
panel_C_2 <- read_fig_proto("PAAD_G12D_FKBP1A_box-plot")

panel_C <- (panel_C_1 + panel_C_2) + plot_layout(heights = c(1, 2))


#### ---- D. PAAD - G12D - FKBP1A genetic graph ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_D <- read_fig_proto("PAAD_G12D_FKBP1A_genetic-graph-plot")



#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (panel_A | panel_B) / (panel_C | panel_D)

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
