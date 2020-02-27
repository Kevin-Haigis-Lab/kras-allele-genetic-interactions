# Build Figure 3

FIGNUM <- 3
VERSION <- 1
FIGFILENAME <- glue("figure_{FIGNUM}_{VERSION}.svg")
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


#### ---- Figure theme ---- ####

#' General theme for Figure 3.
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
# The high-level network plot for the comutation graph for LUAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_LUAD", FIGNUM) +
    theme_graph_fig3() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm")
    ) +
    labs(tag = "a")



#### ---- B. Survival curves of comutated genes ---- ####
# Survival curves of genes comutating with G12C.
# original script: "src/70_15_comutation-survival-analysis.R"

proto_paths <- c("survival_alleleorwt_CHRNB4-G12C-LUAD.rds",
                 "survival_alleleorwt_VN1R2-G12C-LUAD.rds",
                 "survival_alleleorwt_ZNF445-G12C-LUAD.rds",
                 "survival_alleleorwt_ZNF804A-G12C-LUAD.rds")

survival_curves <- purrr::map(
    proto_paths,
    function(x) {
        read_fig_proto(x, FIGNUM) +
        theme_fig3() +
        theme(
            axis.text = element_text(family = "arial", size = 5),
            legend.position = "none"
        ) +
        labs(x = "days")
    })
survival_curves[[1]] <- survival_curves[[1]] + labs(tag = "b")
panel_B <- wrap_plots(survival_curves, ncol = 2)


panel_B_legend <- read_fig_proto("custom_survival_curve_legend", FIGNUM) +
    scale_x_continuous(limits = c(0.7, 4.3)) +
    theme_void() +
    theme(
        legend.position = "none"
    )


#### ---- C. Dot-plot of functional enrichment ---- ####
# A dot plot of the results of functional enrichment in LUAD comutation network.
# original script: "src/20_45_fxnal-enrich-genetic-interactions.R"

panel_C <- read_fig_proto("enrichr_LUAD", FIGNUM) +
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
    theme_fig3() +
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
    labs(
        tag = "c",
        size = expression(-italic("log")[10] ( "adj. p-value" )),
        alpha = "num. of genes"
    )


#### ---- D. Bar-plot of functional enrichment ---- ####
# A bar plot of the frequency of comutation for the enriched functions
# shown in panel B.
# original script: "src/20_46_enriched-functions_bar-plots.R"

panel_D <- read_fig_proto("comut-barplot_LUAD_AllAlleles_AllSources", FIGNUM) +
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
    theme_fig3() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "mm")
    ) +
    labs(
        tag = "d",
        y = expression("" %<-% "reduced | increased" %->% "")
    )



#### ---- E. Labeled MM comutation graph ---- ####
# The comutation graph for MM with every node labeled.
# original script: "src/60_10_MM-specific-oncogenes.R"

panel_E_1 <- read_fig_proto("mm_comut_heatmap", FIGNUM) +
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
panel_E_2 <- read_fig_proto("allele_freq_barplot", FIGNUM) +
    scale_y_continuous(
        breaks = c(10, 50, 200, 700),
        expand = expand_scale(mult = c(0, 0.05)),
        trans = "log10"
    ) +
    theme_fig3() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 5),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, -0.5, 0, "mm")
    ) +
    labs(tag = "e")
panel_E_3 <- read_fig_proto("gene_freq_barplot", FIGNUM) +
    scale_y_continuous(
        breaks = c(50, 100, 200),
        expand = expand_scale(mult = c(0, 0.05)),
    ) +
    theme_fig3() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 5),
        plot.margin = margin(0, 0, 0, -0.5, "mm")
    )

panel_E_design <- "
    11111#
    222223
    222223
    222223
    222223
    222223
"

panel_E <- panel_E_2 + panel_E_1 + panel_E_3 +
    plot_layout(design = panel_E_design)



#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    panel_B <- wrap_elements(
        full = (panel_B) / panel_B_legend + plot_layout(heights = c(10, 1))
    )

    # COMPLETE FIGURE
    top_panels <- (
        (wrap_elements(full = panel_A) / panel_B) |
        (wrap_elements(full = panel_C) / wrap_elements(full = panel_D))
    )

    bottom_panels <- (plot_spacer() | wrap_elements(full = panel_E)) +
        plot_layout(widths = c(1, 1))

    full_figure <- (top_panels | bottom_panels) +
        plot_layout(heights = c(2, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
