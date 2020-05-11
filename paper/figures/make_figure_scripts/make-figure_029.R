# Figure 029. COAD dependency analysis.
# This will be the main figure for the paper and the other cancers will
# be in the supplemental.

FIGNUM <- 29

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 160


theme_fig29 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}



#### ---- A. GSEA ---- ####
# GSEA of the dependencies.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_A <- read_fig_proto("gsea-results-COAD-select.rds") +
    scale_color_gradient2(
        low = synthetic_lethal_pal["down"],
        high = synthetic_lethal_pal["up"],
        guide = guide_colorbar(
            title = "NES",
            title.position = "left",
            label.position = "right",
            barheight = unit(15, "mm"),
            barwidth = unit(2, "mm"),
            order = 1
        )
    ) +
    scale_size_continuous(
        range = c(0.5, 3),
        guide = guide_legend(
            title = "-*log*<sub>10</sub>(adj. p-value)",
            title.position = "left",
            label.position = "right",
            label.hjust = 0,
            keywidth = unit(0, "mm"),
            order = 2
        )
    ) +
    theme_fig29(tag_margin = margin(0, 0, 0, -3, "mm")) +
    theme(
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(-2, 0, 0, 0, "mm"),
        legend.title = element_markdown(angle = 90, vjust = 0.5, hjust = 0.5),
        legend.spacing = unit(0, "mm")
    ) +
    labs(tag = "a")


#### ---- B, C. Ranked heatmaps of GSEA ---- ####
# Two heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

theme_fig29_densityplots <- function(tag_margin_l = -3) {
    theme_classic_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, -2, -1, tag_margin_l, "mm")
        ),
        plot.title = element_text(size = 7, family = "Arial", face = "bold"),
        plot.margin = margin(0.4, 0, -1.5, 0, "mm"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
    )
}

x_label <- expression("" %<-% "greater dep. - ranked by gene effect - less dep." %->% "")

panel_B_density <- read_fig_proto(
        "rankline_COAD_G12V_PID_ERBB4_PATHWAY.svg"
    ) +
    theme_fig29_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        tag = "b",
        title = "ERBB4 pathway"
    )

panel_B <- read_fig_proto(
        "rankplot_COAD_G12V_PID_ERBB4_PATHWAY.svg"
    ) +
    theme_fig29() +
    theme(
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = NA, color = NA)
    ) +
    labs(
        x = x_label
    )


panel_C <- read_fig_proto(
        "rankplot_COAD_G13D_KEGG_OXIDATIVE_PHOSPHORYLATION") +
    theme_fig29() +
    theme(
        plot.title = element_text(size = 7, family = "Arial"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.direction = "none",
        panel.grid = element_blank()
    ) +
    labs(fill = "*KRAS* allele")


panel_C_density <- read_fig_proto(
        "rankline_COAD_G13D_KEGG_OXIDATIVE_PHOSPHORYLATION"
    ) +
    theme_fig29_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        title = "Oxidative phosphorylation",
        tag = "c"
    )

panel_C <- panel_C +
    theme_fig29() +
    theme(
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = NA, color = NA)
    ) +
    labs(
        x = x_label
    )


lbl_alleles <- ggplot_build(panel_C)$plot$data$kras_allele %>%
    unique() %>%
    factor_alleles() %>%
    sort() %>%
    as.character()

panel_BC_legend <- custom_label_legend(
        lbl_alleles,
        gap = 0,
        colors = ifelse(lbl_alleles %in% kras_dark_lbls, "white", "black"),
        y_value = "*KRAS* allele",
        size = 2, fontface = "bold", family = "Arial",
        label.padding = unit(1, "mm"), label.size = unit(0, "mm")
    ) +
    scale_fill_manual(values = short_allele_pal) +
    theme(
        legend.position = "none",
        plot.title = element_blank(),
        axis.text.y = element_markdown(hjust = 0.5, face = "bold",
                                       size = 6, family = "Arial")
    )


#### ---- D. Heatmap of linear model ---- ####
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal in COAD.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_D <- read_fig_proto("COAD_CRISPR_euclidean_ward.D2_pheatmap.rds")
pre_panel_D <- pre_panel_D[[4]]

pre_panel_D_main <- gtable::gtable_filter(pre_panel_D,
                                          "legend",
                                          invert = TRUE)

panel_D <- wrap_elements(plot = pre_panel_D_main) *
    theme_fig29(tag_margin = margin(0, -8, 0, 6, "mm")) %+replace%
    theme(
        plot.margin = margin(-10, -9, -13, -7, "mm")
    ) +
    labs(tag = "d")


panel_D_legend1 <- read_fig_proto(
        "COAD_CRISPR_euclidean_ward.D2_pheatmap_heatpal.rds"
    )

panel_D_legend <- ggplot_build(panel_D_legend1)$plot$data %>%
    mutate(x = 1:n(),
           name = as.character(round(name, 1)),
           name = fct_inorder(name)) %>%
    ggplot(aes(x = name, y = "1")) +
    geom_tile(aes(fill = value), color = NA) +
    scale_fill_identity(guide = FALSE) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_fig34() +
    theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()
    ) +
    labs(title = "scaled dep. score")



#### ---- E. Box plots of select genes ---- ####
# Box-plots of 4 genes manually selected from the heatmap (panel D).
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

panel_E_files <- c(
    "COAD-LIN7C_extra.rds",
    "COAD-TFPT_extra.rds",
    "COAD-STARD9_extra.rds",
    "COAD-KNTC1_extra.rds"
)

panel_E_plots <- as.list(rep(NA, length(panel_E_files)))
names(panel_E_plots) <- panel_E_files
for (f in panel_E_files) {

    if (FALSE) {
        x_axis_text <- element_blank()
    } else {
        x_axis_text <- element_text(size = 5.5, hjust = 0.5, vjust = 1)
    }

    panel_E_plots[[f]] <- read_fig_proto(f) +
        theme_fig29(tag_margin = margin(0, 0, 0, -2, "mm")) %+replace%
        theme(
            plot.title = element_text(size = 6, face = "bold"),
            axis.title.x = element_blank(),
            axis.text.x = x_axis_text,
            legend.position = "none",
            plot.margin = margin(1, 0, 1, 0, "mm")
        ) +
        labs(
            title = str_remove_all(file_sans_ext(f), "COAD-|_extra"),
            y = "dependency score"
        )
}

panel_E <- wrap_plots(panel_E_plots, ncol = 1)
panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")



#### ---- Figure assembly ---- ####

{
    panel_BC_legend_spaced <- (
        plot_spacer() | panel_BC_legend | plot_spacer()
    ) +
        plot_layout(widths = c(1, 5, 1))


    panel_BC <- (panel_B_density / panel_B / panel_C_density / panel_C) +
        plot_layout(heights = c(2, 5, 2, 5))

    column_3 <- (
        panel_E /
        wrap_elements(full = panel_D_legend)
    )

    # COMPLETE FIGURE
    full_figure <- (
        (
            (
                panel_A /
                wrap_elements(full = panel_BC) /
                panel_BC_legend
            ) +
                plot_layout(heights = c(20, 30, 1))
        ) |
        (
            panel_D
        ) |
        (
            panel_E
        )
    ) +
        plot_layout(widths = c(3, 8, 2.5))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
