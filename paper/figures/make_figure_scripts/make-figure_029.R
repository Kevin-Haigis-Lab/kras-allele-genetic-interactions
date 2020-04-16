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

x_label <- expression("" %<-% "greater dep. - ranked gene effect - less dep." %->% "")

panel_B_density <- read_fig_proto(
        "rankline_COAD_G12V_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"
    ) +
    theme_fig29_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        tag = "b",
        title = "Respiratory electron transport"
    )

panel_B <- read_fig_proto(
        "rankplot_COAD_G12V_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT.rds"
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
        "rankplot_COAD_G13D_REACTOME_COMPLEMENT_CASCADE.rds") +
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
        "rankline_COAD_G13D_REACTOME_COMPLEMENT_CASCADE.rds"
    ) +
    theme_fig29_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        title = "Complement cascade",
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


lbl_alleles <- ggplot_build(panel_C)$plot$data$allele %>%
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

pre_panel_D <- read_fig_proto("COAD_CRISPR_manhattan_ward.D2_pheatmap.rds")
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


prep_pheatmap_legend <- function(name) {
    read_fig_proto(name) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_fig29() +
        theme(
            plot.title = element_text(size = 6, family = "Arial"),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 6, family = "Arial"),
            plot.background = element_rect(fill = NA, color = NA),
            panel.background = element_rect(fill = NA, color = NA)
        )
}

prep_pheatmap_colorbar <- function(name) {
    read_fig_proto(name) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_fig29() +
            theme(
                plot.title = element_text(size = 6, family = "Arial"),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 6, family = "Arial"),
                plot.background = element_rect(fill = NA, color = NA),
                panel.background = element_rect(fill = NA, color = NA)
            )
}

panel_D_legend1 <- prep_pheatmap_colorbar(
        "COAD_CRISPR_manhattan_ward.D2_pheatmap_heatpal.rds"
    ) +
    labs(title = "scaled\ndep. score")

panel_D_legend2 <- prep_pheatmap_legend(
        "COAD_CRISPR_manhattan_ward.D2_pheatmap_allelepal.rds"
    ) +
    labs(title = "allele")

panel_D_legend3 <- prep_pheatmap_legend(
        "COAD_CRISPR_manhattan_ward.D2_pheatmap_clusterpal.rds"
    ) +
    labs(title = "cluster")

panel_D_legend <- panel_D_legend1 | panel_D_legend2 | panel_D_legend3


#### ---- E. Box plots of select genes ---- ####
# Box-plots of 4 genes manually selected from the heatmap (panel D).
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

panel_E_files <- c(
    "COAD-KNTC1.rds",
    "COAD-IDH1.rds",
    "COAD-PIP5K1A.rds",
    "COAD-WDR26.rds"
)

panel_E_plots <- as.list(rep(NA, length(panel_E_files)))
names(panel_E_plots) <- panel_E_files
for (f in panel_E_files) {

    if (f != "COAD-WDR26.rds") {
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
            plot.margin = margin(0, 0, 0, 0, "mm")
        ) +
        labs(
            title = str_remove(file_sans_ext(f), "COAD-")
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
            column_3
        )
    ) +
        plot_layout(widths = c(3, 7, 2.1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
