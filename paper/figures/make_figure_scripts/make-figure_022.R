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
        plot.title = element_text(size = 7, face = "bold",
                                  hjust = 0.1, vjust = 5),
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
        legend.position = c(0.05, 0.05),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6)
    ) +
    labs(tag = "a",
         title = "COAD")


#### ---- B. A priori genes of interest comutation network for COAD ---- ####
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"
panel_B <- read_fig_proto(
        "goi_overlap_genetic_interactions_network_COAD_allLists"
    ) +
    scale_size_manual(
        values = c(big = 1.6, small = 1.5),
        guide = FALSE
    ) +
    theme_graph_fig22(plot_margin = margin(0, 0, 0, 0, "mm")) %+replace%
    theme(
        legend.title = element_markdown(size = 6),
        legend.text = element_text(size = 6)
    ) +
    labs(
        tag = "b",
        edge_width = "-*log*(p-value)",
        title = "COAD"
    )


#### ---- C. Dot-plot of functional enrichment  ---- ####
# A dot plot of the results of functional enrichment in COAD comutation network.
# original script: "src/20_45_fxnal-enrich-genetic-interactions.R"

prepare_enrichr_dotplot <- function(plt) {
    plt +
        scale_size_area(
            max_size = 4,
            guide = guide_legend(
                title.position = "left",
                title.hjust = 0.5,
                order = 10,
                label.vjust = 0,
                label.position = "top"
            )
        ) +
        scale_alpha_continuous(
            range = c(0.3, 1),
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
            axis.text.x = element_text(size = 6),
            legend.title = element_markdown(),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.margin = margin(-1, 3, -1, 3, "mm"),
            legend.box.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            strip.text = element_text(size = 7, hjust = 0.5, face = "bold")
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
        top_bar_limits = NULL,
        top_bar_breaks = integer_breaks(rm_vals = c(0)),
        top_bar_labels = waiver(),
        right_bar_limits = NULL,
        right_bar_breaks = integer_breaks(rm_vals = c(0)),
        right_bar_labels = waiver(),
        tag_margin = margin(0, 0, 0, 0, "mm")
    ) {
    # Top bar plot
    pw[[1]] <- pw[[1]] +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.01)),
            limits = top_bar_limits,
            breaks = top_bar_breaks,
            labels = top_bar_labels
        ) +
        theme(
            axis.text.y = element_text(size = 6, hjust = 1),
            plot.tag = element_text(size = 7,
                                    face = "bold",
                                    margin = tag_margin),
            axis.ticks = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "mm")
        )


    # Middle main tile plot
    pw[[2]] <- pw[[2]] +
        scale_fill_manual(
            values = mod_variant_pal,
            guide = guide_legend(title = "mutation type",
                                 nrow = 2,
                                 label.theme = element_text(size = 5.5),
                                 title.position = "top")
        ) +
        theme(
            axis.text.y = element_text(size = 6, hjust = 1),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            legend.background = element_rect(fill = NULL, color = NULL),
            legend.margin = margin(-5, 0, -5, 0, "mm"),
            legend.spacing.x = unit(1, "mm")
        )

    # Right bar plot
    pw[[3]] <- pw[[3]] +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.08)),
            limits = right_bar_limits,
            breaks = right_bar_breaks,
            labels = right_bar_labels
        ) +
        theme(
            axis.text.x = element_text(size = 6, vjust = 1),
            axis.ticks = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "mm")
        )


    oncoplot_design <- "
        11111#
        11111#
        222223
        222223
        222223
        222223
        222223
    "

    pw <- pw + plot_layout(design = oncoplot_design)

    return(pw)
}


# Set the legend position for the oncoplots to "none".
remove_oncoplot_legend <- function(pw) {
    # Only panel 2 has a legend.
    pw[[2]] <- pw[[2]] + theme(legend.position = "none")
    return(pw)
}


panel_D <- read_fig_proto("COAD_G12D_comutation_oncostrip_select_NO-KRAS-ANNO")
panel_D <- adjust_oncoplot_theme(panel_D,
                                 right_bar_limits = c(0, 510),
                                 right_bar_breaks = c(100, 200, 300, 400),
                                 right_bar_labels = c("", "200", "", "400"),
                                 tag_margin = margin(-3, 0, 0, -3.5, "mm"))
panel_D[[1]] <- panel_D[[1]] + labs(tag = "d")

# Legend for panels d-g
# panel_D_leg_1 <- ggpubr::as_ggplot(cowplot::get_legend(panel_D[[2]]))
# panel_D_leg <- (plot_spacer() / panel_D_leg_1 / plot_spacer()) +
#     plot_layout(heights = c(2, 1, 2))
# panel_D_leg <- wrap_elements(full = panel_D_leg)
panel_D <- remove_oncoplot_legend(panel_D)


#### ---- E. Oncoplot ---- ####
# A rainfall plot of select reduced comutation interactions with G12D in COAD.
# original script: "src/20_50_rainfall-plots.R"

panel_E <- read_fig_proto("COAD_G12D_exclusivity_oncostrip_select_NO-KRAS-ANNO")
panel_E <- adjust_oncoplot_theme(panel_E,
                                 top_bar_breaks = c(1:4),
                                 top_bar_labels = c("", "2", "", "4"),
                                 right_bar_limits = c(0, 510),
                                 right_bar_breaks = c(100, 200, 300, 400),
                                 right_bar_labels = c("", "200", "", "400"),
                                 tag_margin = margin(0, 0, 0, -3.5, "mm"))
panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")
panel_E <- remove_oncoplot_legend(panel_E)


#### ---- Oncoplot mutation legend ---- ####

extract_variants_from_oncoplots <- function(p) {
    ggplot_build(p[[2]])$plot$data$Variant_Classification
}

calc_starts <- function(x, gap = 0, starting_pt = 0) {
    starts <- accumulate(x, .init = starting_pt, ~ .x + .y + gap)
    starts[-length(starts)]
}

variants <- unique(c(extract_variants_from_oncoplots(panel_D),
                     extract_variants_from_oncoplots(panel_E)))

panel_D_leg_df <- tibble(lbl = variants) %>%
    mutate(
        len = str_length(lbl),
        start = calc_starts(len, gap = 0),
        end = len + start,
        mid = (start + end) / 2
    )

panel_D_leg <- panel_D_leg_df %>%
    ggplot(aes(x = mid, y = 1)) +
    geom_text(aes(label = lbl, color = lbl),
              family = "Arial", size = 1.8, fontface = "bold") +
    scale_color_manual(values = mod_variant_pal, guide = FALSE) +
    scale_x_continuous(
        limits = c(min(panel_D_leg_df$start), max(panel_D_leg_df$end)),
        expand = c(0, 0)
    ) +
    theme_void() +
    theme(
        plot.title = element_text(size = 6, hjust = 0),
        axis.title = element_blank(),
        axis.text = element_blank()
    ) +
    labs(title = "mutation type")


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
        plot.title = element_text(hjust = 0, face = "bold",
                                  size = 7, vjust = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_textbox_simple(
            size = 6,
            hjust = 0.5,
            vjust = 0,
            padding = margin(0, 0, 0, 0),
            margin = margin(0, 0, -10, 0),
            halign = 0.5,
            orientation = "left-rotated",
        ),
        axis.text.x = element_blank(),
        plot.margin = margin(6, 0, -0.5, 1, "mm")
    ) +
    labs(tag = "f",
         y = "num. tumor<br>samples<br>(*log*<sub>10</sub>)",
         title = "MM")
panel_F_3 <- read_fig_proto("gene_freq_barplot_TRUNCATED") +
    scale_y_continuous(
        breaks = c(25, 100, 200),
        expand = expansion(mult = c(0, 0.05)),
    ) +
    theme_fig22() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 6),
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

    panels_DE <- (panel_D / panel_E / wrap_elements(full = panel_D_leg)) +
        plot_layout(heights = c(15, 15, 1)) +
        plot_annotation(
            title = "COAD",
            theme = theme(plot.title = element_text(size = 7,
                                                    face = "bold",
                                                    hjust = 0.1,
                                                    vjust = -1))
        )
    panels_DE <- wrap_elements(full = panels_DE)

    row_1 <- wrap_elements(full = panel_A | panel_B)

    # COMPLETE FIGURE
    full_figure <- (
        row_1 /
        panel_C /
        ((panels_DE | panel_F) + plot_layout(widths = c(2, 3)))
    ) +
        plot_layout(heights = c(4, 3, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
