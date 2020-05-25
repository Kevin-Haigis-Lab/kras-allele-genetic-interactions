# Figure 046. #> Second half of results from explaining allele-specific
# genetic dependency with comutation interactions.

FIGNUM <- 46

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig46 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}



theme_coefplot <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_fig46(tag_margin = tag_margin) %+replace%
    theme(
        plot.title = element_markdown(size = 7,
                                      hjust = 0,
                                      vjust = 0.5,
                                      valign = 1,
                                      face = "bold",
                                      family = "Arial"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown()
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



#### ---- A. COAD - G12D - SRSF5 ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("COAD_G12D_SRSF5_coef-plot") +
    scale_y_discrete(expand = expansion(mult = c(0.3, 0.3))) +
    theme_coefplot() +
    labs(tag = "a",
         title = "COAD")
panel_A_2 <- read_fig_proto("COAD_G12D_SRSF5_box-plot") +
    theme_fig46() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        legend.title = element_markdown(),
        legend.spacing.y = unit(1, "mm"),
        legend.key.height = unit(3, "mm")
    )
panel_A_3 <- read_fig_proto("COAD_G12D_SRSF5_comut-graph-plot")


panel_A <- (panel_A_1 / panel_A_2 / panel_A_3) +
    plot_layout(heights = c(3, 9, 1))



#### ---- B. PAAD - G12R - KIAA1257 ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12R_KIAA1257_coef-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    theme_coefplot() +
    labs(tag = "b",
         title = "PAAD")
panel_B_2 <- read_fig_proto("PAAD_G12R_KIAA1257_box-plot")+
    theme_fig46() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        legend.position = "none"
    )
panel_B_3 <- read_fig_proto("PAAD_G12R_KIAA1257_G12R-DNAH5_line-plot") +
    theme_fig46_void() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank()
    )


panel_B_layout <- c(
    area(t = 1, l = 1, b = 4, r = 10),
    area(t = 5, l = 1, b = 20, r = 10),
    area(t = 6.5, l = 5.5, b = 8, r = 10)
)

panel_B <- panel_B_1 + panel_B_2 + panel_B_3 +
    plot_layout(design = panel_B_layout)


#### ---- C. PAAD - G12D - FKBP1A ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_C_1 <- read_fig_proto("PAAD_G12D_FKBP1A_coef-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.2, 0.05))) +
    theme_coefplot() +
    labs(tag = "c",
         title = "PAAD")
panel_C_2 <- read_fig_proto("PAAD_G12D_FKBP1A_box-plot") +
theme_fig46() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        legend.title = element_markdown(),
        legend.position = "bottom",
        legend.spacing.y = unit(1, "mm"),
        legend.key.height = unit(3, "mm")
    )

panel_C <- (panel_C_1 + panel_C_2) + plot_layout(heights = c(2, 7))


#### ---- D. PAAD - G12D - FKBP1A genetic graph ---- ####
# original script: "src/40_63_explore-synlet-comut.R"

panel_D <- read_fig_proto("PAAD_G12D_FKBP1A_genetic-graph-plot") +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme_fig46_void() +
    theme(
        plot.margin = margin(3, 3, 3, 3, "mm"),
        legend.position = "bottom",
        legend.title = element_markdown(),
        legend.spacing.y = unit(1, "mm"),
        legend.key.height = unit(3, "mm")
    ) +
    labs(tag = "d")



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
