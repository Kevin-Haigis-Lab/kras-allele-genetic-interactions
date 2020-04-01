# Figure 003. Main comutation figure for LUAD and MM.

FIGNUM <- 3

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


#### ---- Figure theme ---- ####

#' General theme for Figure 3.
theme_fig3 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
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
# The high-level network plot for the comutation graph for LUAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_LUAD") +
    theme_graph_fig3() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = c(0.15, 0.9)
    ) +
    labs(tag = "a")



#### ---- B. A priori genes of interest comutation network for LUAD ---- ####
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"
panel_B <- read_fig_proto(
        "goi_overlap_genetic_interactions_network_LUAD_kegg"
    ) *
    theme_graph_fig3() %+replace%
    theme(
        legend.title = element_markdown()
    ) +
    labs(
        tag = "b",
        edge_width = "-*log*(p-value)"
    )


#### ---- C. Dot-plot of functional enrichment ---- ####
# A dot plot of the results of functional enrichment in LUAD comutation network.
# original script: "src/20_45_fxnal-enrich-genetic-interactions.R"

panel_C <- read_fig_proto("enrichr_LUAD") +
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
    theme_fig3(margin(-1, 0, 0, 0, "mm")) +
    theme(
        plot.title = element_blank(),
        axis.title = element_blank(),
        legend.title = element_markdown(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(-1, 0, -1, 0, "mm"),
        legend.box.background = element_rect(fill = NA, color = NA)
    ) +
    labs(
        tag = "c",
        size = "-*log*<sub>10</sub>(adj. p-value)",
        alpha = "num. of genes"
    )


#### ---- D. Bar-plot of functional enrichment ---- ####
# A bar plot of the frequency of comutation for the enriched functions
# shown in panel B.
# original script: "src/20_46_enriched-functions_bar-plots.R"

panel_D <- read_fig_proto("comut-barplot_LUAD_AllAlleles_AllSources") +
    scale_fill_manual(
        values = comut_updown_pal,
        guide = guide_legend(
            title.position = "top"
        )
    ) +
    scale_alpha_manual(
        values = c("other" = 0.4, "allele" = 0.95),
        guide = guide_legend(
            title.position = "top"
        )
    ) +
    scale_y_continuous(
        label = abs
    ) +
    theme_fig3() %+replace%
    theme(
        legend.title = element_markdown(),
        legend.spacing.x = unit(1, "mm"),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "mm"),
        axis.text = element_text(size = 5, family = "arial")
    ) +
    labs(
        tag = "d",
        y = expression("" %<-% "reduced | increased" %->% ""),
        alpha = "*KRAS*"
    )


#### ---- A LINE ---- ####

SEPARATING_LINE <- tibble(x = c(1, 2), y = c(1, 1)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_path(color = "grey70", size = 0.3, lineend = "round") +
    scale_x_continuous(limits = c(1, 2)) +
    theme_void() +
    theme(
        plot.margin = margin(0, 0, 0, 0, "mm")
    )


#### ---- E. High-level comutation network for MM (labeled) ---- ####
# The labeled, high-level network plot for the comutation graph for MM
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_E <- read_fig_proto("genetic_interaction_network_labeled_MM") +
    # scale_color_manual(values = short_allele_pal, guide = FALSE) +
    theme_graph_fig3() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = c(0.1, 0.1)
    ) +
    labs(tag = "e")


#### ---- F. Labeled MM comutation graph ---- ####
# The comutation graph for MM with every node labeled.
# original script: "src/60_10_MM-specific-oncogenes.R"

panel_F_1 <- read_fig_proto("mm_comut_heatmap") +
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "tomato",
        mid = "grey90",
        midpoint = 0.10,
        guide = guide_colorbar(
            barheight = unit(2, "mm")
        )
    ) +
    theme_fig3() +
    theme(
        legend.position = "none",
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
panel_F_2 <- read_fig_proto("allele_freq_barplot") +
    scale_y_continuous(
        breaks = c(10, 50, 200, 700),
        expand = expansion(mult = c(0, 0.05)),
        trans = "log10"
    ) +
    theme_fig3(margin(-1.9, 0, 0, 0, "mm")) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_textbox_simple(
            size = 5,
            hjust = 0.5,
            vjust = 0,
            padding = margin(0, 0, 0, 0),
            margin = margin(0, 0, -10, 0),
            halign = 0.5,
            orientation = "left-rotated",
        ),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, -0.5, 0, "mm")
    ) +
    labs(tag = "f",
         y = "num. tumor<br>samples<br>(*log*<sub>10</sub>)")
panel_F_3 <- read_fig_proto("gene_freq_barplot") +
    scale_y_continuous(
        breaks = c(25, 100, 200),
        expand = expansion(mult = c(0, 0.05)),
    ) +
    theme_fig3() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 5),
        plot.margin = margin(0, 0, 0, -0.5, "mm")
    )

panel_F_design <- "
    11111#
    222223
    222223
    222223
    222223
    222223
"

panel_F <- panel_F_2 + panel_F_1 + panel_F_3 +
    plot_layout(design = panel_F_design)



#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    # COMPLETE FIGURE
    top_panels <- (
        (wrap_elements(full = panel_A) / panel_B) |
        (wrap_elements(full = panel_C) / wrap_elements(full = panel_D))
    )

    bottom_panels <- (
            wrap_elements(full = panel_E) | wrap_elements(full = panel_F)
        ) +
            plot_layout(widths = c(1, 1))

    full_figure <- (top_panels / SEPARATING_LINE / bottom_panels) +
        plot_layout(heights = c(20, 1, 10))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
