# Figure 002. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 2

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


#### ---- Figure theme ---- ####

#' General theme for Figure 2.
theme_fig2 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm"))
    )
}


#' Special theme for graphs from 'ggraph'.
theme_graph_fig2 <- function() {
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


#' Adjust the oncoplots.
adjust_oncoplot_theme <- function(
        pw,
        right_bar_limits = NULL,
        right_bar_breaks = integer_breaks(rm_vals = c(0)),
        tag_margin = margin(0, 0, 0, 0, "mm")
    ) {
    # Top bar plot
    pw[[1]] <- pw[[1]] +
        theme(
            axis.text.y = element_text(size = 5, hjust = 1),
            plot.tag = element_text(size = 7,
                                    face = "bold",
                                    margin = tag_margin)
        )


    # Middle main tile plot
    pw[[2]] <- pw[[2]] +
        scale_fill_manual(
            values = mod_variant_pal,
            guide = guide_legend(title = "mutation type",
                                 ncol = 1,
                                 title.position = "top")
        ) +
        theme(
            axis.text.y = element_text(size = 6, hjust = 1)
        )

    # Right bar plot
    pw[[3]] <- pw[[3]] +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.08)),
            limits = right_bar_limits,
            breaks = right_bar_breaks
        ) +
        theme(
            axis.text.x = element_text(size = 5, vjust = 1)
        )


    # Clinical feature bar (KRAS allele)
    pw[[4]] <- pw[[4]] +
        scale_fill_manual(
            values = short_allele_pal,
            guide = guide_legend(
                title = "*KRAS*<br>allele",
                ncol = 1,
                label.hjust = 0,
                label.position = "right",
            )
        ) +
        theme(
            legend.spacing.x = unit(2, "mm"),
            legend.title = element_markdown()
        )

    return(pw)
}


#' Set the legend position for the oncoplots to "none".
remove_oncoplot_legend <- function(pw) {
    # Only panels 2 and 4 have legends.
    for (i in c(2, 4)) {
        pw[[i]] <- pw[[i]] + theme(legend.position = "none")
    }
    return(pw)
}




#### ---- A. High-level comutation network for COAD ---- ####
# The high-level network plot for the comutation graph for COAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"
panel_A <- read_fig_proto("genetic_interaction_network_COAD") +
    scale_color_manual(values = short_allele_pal, guide = FALSE) +
    theme_graph_fig2() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm"),
        legend.position = c(0.1, 0.1)
    ) +
    labs(tag = "a")


#### ---- B. A priori genes of interest comutation network for COAD ---- ####
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"
panel_B <- read_fig_proto(
        "goi_overlap_genetic_interactions_network_COAD_allLists"
    ) *
    theme_graph_fig2() %+replace%
    theme(
        legend.title = element_markdown()
    ) +
    labs(
        tag = "b",
        edge_width = "-*log*(p-value)"
    )


#### ---- C. Dot-plot of functional enrichment  ---- ####
# A dot plot of the results of functional enrichment in COAD comutation network.
# original script: "src/20_45_fxnal-enrich-genetic-interactions.R"
panel_C <- read_fig_proto("enrichr_COAD") +
    scale_size_continuous(
        range = c(0, 3),
        guide = guide_legend(
            title.position = "left",
            title.hjust = 0.5,
            order = 10,
            label.vjust = -11,
            label.position = "top"
        )
    ) +
    scale_alpha_continuous(
        range = c(0, 1),
        guide = guide_legend(
            title.position = "left",
            title.hjust = 0,
            order = 20,
            label.vjust = -11,
            label.position = "top"
        )
    ) +
    theme_fig2() +
    theme(
        plot.title = element_blank(),
        axis.title = element_blank(),
        legend.title = element_markdown(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(-1, 1, -1, -2, "mm"),
        legend.box.background = element_rect(fill = NA, color = NA)
    ) +
    labs(
        tag = "c",
        size = "-*log*<sub>10</sub>(adj. p-value)",
        alpha = "num. of genes"
    )


#### ---- D. Oncoplot ---- ####
# A rainfall plot of select increased comutation interactions with G12D in COAD.
# original script: "src/20_50_rainfall-plots.R"

panel_D <- read_fig_proto("COAD_G12D_comutation_oncostrip_select")
panel_D <- adjust_oncoplot_theme(panel_D, c(0, 450), c(50, 200, 400),
                                 margin(0, 0, 0, -3.5, "mm"))
panel_D[[1]] <- panel_D[[1]] + labs(tag = "d")

# Legend for panels d-g
panel_D_leg_1 <- ggpubr::as_ggplot(cowplot::get_legend(panel_D[[2]]))
panel_D_leg_2 <- ggpubr::as_ggplot(cowplot::get_legend(panel_D[[4]]))
panel_D_leg <- (
        plot_spacer() / panel_D_leg_1 /
        plot_spacer() / panel_D_leg_2 /
        plot_spacer()
    ) +
    plot_layout(heights = c(20, 10, 5, 10, 20))
panel_D_leg <- wrap_elements(full = panel_D_leg)
panel_D <- remove_oncoplot_legend(panel_D)


#### ---- E. Oncoplot ---- ####
# A rainfall plot of select reduced comutation interactions with G12D in COAD.
# original script: "src/20_50_rainfall-plots.R"

panel_E <- read_fig_proto("COAD_G12D_exclusivity_oncostrip_select")
panel_E <- adjust_oncoplot_theme(panel_E, c(0, 450), c(50, 200, 400))
panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")
panel_E <- remove_oncoplot_legend(panel_E)


#### ---- F. Oncoplot ---- ####
# A rainfall plot of select increased comutation interactions with G12D in COAD.
# original script: "src/20_50_rainfall-plots.R"

panel_F <- read_fig_proto("COAD_G12V_comutation_oncostrip_select")
panel_F <- adjust_oncoplot_theme(panel_F, c(0, 1150), c(100, 500, 1000),
                                 margin(0, 0, 0, -3.5, "mm"))
panel_F[[1]] <- panel_F[[1]] + labs(tag = "f")
panel_F <- remove_oncoplot_legend(panel_F)


#### ---- G. Oncoplot ---- ####
# A rainfall plot of select reduced comutation interactions with G12D in COAD.
# original script: "src/20_50_rainfall-plots.R"

panel_G <- read_fig_proto("COAD_G12V_exclusivity_oncostrip_select")
panel_G <- adjust_oncoplot_theme(panel_G, c(0, 490), c(50, 200, 400))
panel_G[[1]] <- panel_G[[1]] + labs(tag = "g")
panel_G <- remove_oncoplot_legend(panel_G)


#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    row_1 <- (wrap_elements(full = panel_A) | panel_B | panel_C) +
        plot_layout(widths = c(2, 2, 1))

    rows_2_3 <- (
            (
                (panel_D | panel_E) / (panel_F | panel_G )
            ) |
            panel_D_leg
        ) +
        plot_layout(widths = c(20, 3))

    # COMPLETE FIGURE
    full_figure <- (row_1) / (rows_2_3) +
        plot_layout(heights = c(2, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
