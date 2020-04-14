# Figure 026. #> Supplemental comutation information for PAAD

FIGNUM <- 26

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")
FIG_DIMENSIONS$height <- 90

theme_fig26 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


theme_graph_fig26 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm")),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "bottom",
        legend.title = element_text(size = 7, hjust = 0.5),
        legend.text = element_text(size = 6, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}


#### ---- A. Network of calcium signaling ---- ####
# A network plot of the calcium signaling pathway enriched in the G12V
# comutation network.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_A <- read_fig_proto(
        "PAAD_G12V_GO-Biological-Process-2018_calcium-ion-transport.rds"
    ) +
    scale_color_manual(
        values = c(comut_updown_pal,
                   "none" = "grey70",
                   "in_geneset" = "grey40"),
        guide = guide_legend(
            title = "comutation",
            title.position = "top",
            title.hjust = 0.0,
            title.theme = element_text(size = 6, face = "bold"),
            label.position = "right",
            label.hjust = 0.0,
            legend.text = element_text(size = 6, hjust = 0),
            nrow = 1,
            keywidth = unit(1, "mm"),
            keyheight = unit(1, "mm")
        )
    ) +
    theme_graph_fig26() %+replace%
    theme(
        legend.spacing = unit(0, "mm"),
        legend.position = c(0.2, 0.00),
        legend.margin = margin(-2, 0, 2, 0, "mm"),
    )  +
    labs(tag = "a")


#### ---- B. Distribution of comutation events ---- ####
# A dot-plot of some selected enriched functions from the comutation network.
# original script: "src/20_48_enriched-functions_compare-functions_heatmaps.R"

panel_B <- read_fig_proto("comparison-heatmap_PAAD-1.rds") &
    theme_fig26()

panel_B[[1]] <- panel_B[[1]] +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, hjust = 1.0),
        legend.position = "none"
    ) +
    labs(tag = "b")

panel_B[[2]] <- panel_B[[2]] +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(4, "mm"),
        legend.title = element_markdown(size = 6),
        legend.text = element_text(size = 6),
        axis.ticks.y = element_line(size = 0.1)
    ) +
    labs(y = "distribution of comutation events",
         fill = "*KRAS* allele")



#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- panel_A | panel_B +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
