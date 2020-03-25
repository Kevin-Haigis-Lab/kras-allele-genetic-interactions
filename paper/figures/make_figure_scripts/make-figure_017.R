# Figure 017. The main comutation figure for PAAD.

FIGNUM <- 17

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig17 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


#' Special theme for graphs from 'ggraph'.
theme_graph_figS12 <- function() {
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


#### ---- A. High-level comutation network ---- ####
# A The high-level overview of the comutation network for the *KRAS* alleles
# in PAAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_PAAD") +
    scale_edge_color_manual(
        values = comut_updown_pal,
        guide = guide_legend(
            title = "comutation",
            title.position = "top",
            keywidth = unit(2, "mm"),
            keyheight = unit(1, "mm"),
            nrow = 2,
            label.position = "top",
            order = 1
        )
    ) +
    scale_color_manual(
        values = short_allele_pal,
        na.value = NA,
        guide = FALSE
    ) +
    theme_graph_figS12() +
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = c(0.05, 0.1),
        legend.box = "horizontal",
        legend.margin = margin(-6, 0, 0, 0, "mm")
    )

panel_A <- wrap_elements(full = panel_A) +
    labs(tag = "a") +
    theme_fig17()


#### ---- B. Labeled comutation network of genes in a priori lists ---- ####
# A The labeled comutation network for the *KRAS* alleles only including genes
# in a priori selected gene sets.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"

panel_B <- read_fig_proto(
        "goi_overlap_genetic_interactions_network_PAAD_allLists",
    ) +
    theme_graph_figS12() +
    theme(
        legend.position = "bottom",
        legend.margin = margin(-2, 0, 2, 0, "mm"),
        legend.title = element_markdown()
    ) +
    labs(edge_width = "-*log*(p-value)")

panel_B <- wrap_elements(full = panel_B) +
    labs(tag = "b") +
    theme_fig17()



#### ---- C. Bar plots of log(OR) of select genes ---- ####
# Some genes have two different comutation interactions with multiple KRAS
# alleles. These are bar plots of the log(OR) of a selection of those genes.
# original script: "src/20_41_disagreeing-interactions_logOR-barplot.R"

panel_C <- read_fig_proto(
        "log-odds-ratio_barplot_PAAD.rds",
    ) +
    facet_wrap(~ hugo_symbol, scales = "free", ncol = 1) +
    theme_fig17() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text = element_text(face = "bold")
    ) +
    labs(tag = "c", y = "*log*(OR)")



#### ---- D. Dot-plot of enriched functions ---- ####
# A dot-plot of some selected enriched functions from the comutation network.
# original script: "src/20_41_disagreeing-interactions_logOR-barplot.R"

panel_D <- read_fig_proto("enrichr_PAAD.rds") +
    scale_size_continuous(
        range = c(0, 6),
        guide = guide_legend(title.position = "top",
                             title.hjust = 0.5,
                             title.vjust = 0.0,
                             label.position = "top",
                             label.vjust = -13,
                             nrow = 1,
                             keywidth = unit(1.5, "mm"),
                             keyheight = unit(0, "mm"),
                             order = 10)
    ) +
    scale_alpha_continuous(
        range = c(0.1, 1),
        breaks = c(2, 4, 6, 8),
        guide = guide_legend(title.position = "top",
                             title.hjust = 0.5,
                             label.position = "top",
                             nrow = 1,
                             keywidth = unit(2, "mm"),
                             keyheight = unit(2, "mm"),
                             order = 20)
    ) +
    theme_fig17() +
    theme(
        legend.margin = margin(-6, 0, -3, 0, "mm"),
        legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.box = "horizontal",
        legend.title = element_markdown(),
    ) +
    labs(
        tag = "d",
        size = "-*log*<sub>10</sub>(adj. p-value)",
        alpha = "num. of genes"
    )


#### ---- E. Network of calcium signaling ---- ####
# A network plot of the calcium signaling pathway enriched in the G12V
# comutation network.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_E <- read_fig_proto(
        "PAAD_G12V_GO-Biological-Process-2018_calcium-ion-transport.rds",
    ) +
    scale_color_manual(
        values = c(comut_updown_pal,
                   "none" = "grey70",
                   "in_geneset" = "grey40"),
        guide = guide_legend(
            name = "comutation",
            title.position = "top",
            title.hjust = 0.0,
            label.position = "right",
            label.hjust = 0.0,
            nrow = 1,
            keywidth = unit(1, "mm"),
            keyheight = unit(1, "mm")
        )
    ) +
    theme_graph_figS12() %+replace%
    theme(
        legend.spacing = unit(0, "mm"),
        legend.position = c(0.87, 0.01),
        legend.margin = margin(-2, 0, 2, 0, "mm"),
        legend.text = element_text(size = 5, hjust = 0)
    )

panel_E <- wrap_elements(full = panel_E) +
    labs(tag = "e") +
    theme_fig17()


#### ---- F. Distribution of comutation events ---- ####
# A dot-plot of some selected enriched functions from the comutation network.
# original script: "src/20_48_enriched-functions_compare-functions_heatmaps.R"

panel_F <- read_fig_proto("comparison-heatmap_PAAD-1.rds") &
    theme_fig17()

panel_F[[1]] <- panel_F[[1]] +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, hjust = 1.0),
        legend.position = "none"
    ) +
    labs(tag = "f")

panel_F[[2]] <- panel_F[[2]] +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 5),
        legend.key.size = unit(2, "mm"),
        legend.title = element_blank(),
        axis.ticks.y = element_line(size = 0.1)
    )


#### ---- Figure assembly ---- ####

{
    set.seed(1)

    # COMPLETE FIGURE
    full_figure <-
        wrap_elements(
            full = (
                panel_A + panel_B + panel_C + plot_layout(widths = c(5, 5, 1))
            )
        ) /
        wrap_elements(
            full = (
                panel_D + (
                    (panel_E) / panel_F
                ) +
                plot_layout(widths = c(1, 2))
            )
        ) +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
