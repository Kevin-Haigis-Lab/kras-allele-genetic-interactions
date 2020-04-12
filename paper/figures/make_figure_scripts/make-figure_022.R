# Figure 022. Main comutation interaction figure.

FIGNUM <- 22

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig22 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


# Special theme for graphs from 'ggraph'.
theme_graph_fig22 <- function(plot_margin = margin(0, 0, 0, 0, "mm")) {
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
        plot.margin = plot_margin
    )
}



#### ---- A. High-level comutation network for COAD ---- ####
# The high-level network plot for the comutation graph for COAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"
panel_A <- read_fig_proto("genetic_interaction_network_COAD") +
    theme_graph_fig22() %+replace%
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
    theme_graph_fig22(plot_margin = margin(0, 0, 0, 0, "mm")) %+replace%
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

prepare_enrichr_dotplot <- function(plt) {
    plt +
        scale_size_continuous(
            range = c(0, 4),
            guide = guide_legend(
                title.position = "left",
                title.hjust = 0.5,
                order = 10,
                label.vjust = 0,
                label.position = "top"
            )
        ) +
        scale_alpha_continuous(
            range = c(0, 1),
            guide = guide_legend(
                title.position = "left",
                title.hjust = 0,
                order = 20,
                label.vjust = 0,
                label.position = "top"
            )
        ) +
        theme_fig22() +
        theme(
            plot.title = element_blank(),
            axis.title = element_blank(),
            legend.title = element_markdown(),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.margin = margin(-1, 3, -1, 3, "mm"),
            legend.box.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(0, 0, 0, 0, "mm")
        ) +
        labs(
            tag = "c",
            alpha = "-*log*<sub>10</sub>(adj. p-value)",
            size = "num. of genes"
        )
}

panel_C <- read_fig_proto("enrichr_all-cancers-faceted") %>%
    prepare_enrichr_dotplot()


#### ---- D. Oncoplot ---- ####
# A rainfall plot of select increased comutation interactions with G12D in COAD.
# original script: "src/20_50_rainfall-plots.R"

# Adjust the oncoplots.
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


# Set the legend position for the oncoplots to "none".
remove_oncoplot_legend <- function(pw) {
    # Only panels 2 and 4 have legends.
    for (i in c(2, 4)) {
        pw[[i]] <- pw[[i]] + theme(legend.position = "none")
    }
    return(pw)
}


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


#### ---- F. Labeled MM comutation graph ---- ####
# The comutation graph for MM with every node labeled.
# original script: "src/60_10_MM-specific-oncogenes.R"

panel_F_1 <- read_fig_proto("mm_comut_heatmap_TRUNCATED") +
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "tomato",
        mid = "grey90",
        midpoint = 0.10,
        guide = guide_colorbar(
            barheight = unit(2, "mm")
        )
    ) +
    theme_fig22() +
    theme(
        legend.position = "none",
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 1, "mm")
    )
panel_F_2 <- read_fig_proto("allele_freq_barplot_TRUNCATED") +
    scale_y_continuous(
        breaks = c(10, 50, 200, 700),
        expand = expansion(mult = c(0, 0.05)),
        trans = "log10"
    ) +
    theme_fig22(margin(-1.9, 0, 0, 0, "mm")) +
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
        plot.margin = margin(3, 0, -0.5, 1, "mm")
    ) +
    labs(tag = "f",
         y = "num. tumor<br>samples<br>(*log*<sub>10</sub>)")
panel_F_3 <- read_fig_proto("gene_freq_barplot_TRUNCATED") +
    scale_y_continuous(
        breaks = c(25, 100, 200),
        expand = expansion(mult = c(0, 0.05)),
    ) +
    theme_fig22() +
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
"

panel_F <- panel_F_2 + panel_F_1 + panel_F_3 +
    plot_layout(design = panel_F_design)



#### ---- Figure assembly ---- ####

{
    set.seed(0)

    row_1 <- wrap_elements(full = panel_A | panel_B)

    panels_DE <- wrap_elements(full = panel_D / panel_E)
    # row_3 <- wrap_elements(full = panels_DE | panel_F)

    # COMPLETE FIGURE
    full_figure <- row_1 / panel_C / (panels_DE | panel_F) +
        plot_layout(heights = c(3, 2, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
