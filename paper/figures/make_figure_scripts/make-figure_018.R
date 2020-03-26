# Figure 018. DepMap analysis of PAAD cell lines.

FIGNUM <- 18

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig18 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. GSEA dot-plot ---- ####
# A dot-plot showing some interesting gene sets enriched for the KRAS alleles.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_A <- read_fig_proto("gsea-results-PAAD-select") +
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
    theme_fig18(tag_margin = margin(0, 0, 0, -3, "mm")) +
    theme(
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(-2, 0, 0, 0, "mm"),
        legend.title = element_markdown(angle = 90, vjust = 0.5, hjust = 0.5),
        legend.spacing = unit(0, "mm")
    ) +
    labs(tag = "a")


#### ---- B. GSEA ranked-heatmap (1) ---- ####
# A heatmap of an enriched gene set with the cell lines ranked by their
# dependency score for each gene. A density plot along the top helps
# highlight the trend.
# original script: "src/10_37_gsea-depmap-analysis.R"

theme_fig18_densityplots <- function(tag_margin_l = -3) {
    theme_classic_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, -2, -1, tag_margin_l, "mm")
        ),
        plot.title = element_text(size = 6, family = "Arial"),
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
        "rankline_PAAD_G12D_REACTOME_G_ALPHA_12_13_SIGNALLING_EVENTS"
    ) +
    theme_fig18_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        tag = "b",
        title = "G alpha (12/13) signalling events"
    )

panel_B <- read_fig_proto(
        "rankplot_PAAD_G12D_REACTOME_G_ALPHA_12_13_SIGNALLING_EVENTS"
    ) +
    theme_fig18() +
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



#### ---- C. GSEA ranked-heatmap (2) ---- ####
# A heatmap of an enriched gene set with the cell lines ranked by their
# dependency score for each gene. A density plot along the top helps
# highlight the trend.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_C <- panel_C_density <- read_fig_proto(
        "rankline_PAAD_G12R_REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT"
    ) +
    theme_fig18_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        tag = "c",
        title = "G2/M DNA damage checkpoint"
    )

panel_C <- read_fig_proto(
        "rankplot_PAAD_G12R_REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT"
    ) +
    theme_fig18() +
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


#### ---- D. GSEA ranked-heatmap (3) ---- ####
# A heatmap of an enriched gene set with the cell lines ranked by their
# dependency score for each gene. A density plot along the top helps
# highlight the trend.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_D_density <- read_fig_proto(
        "rankline_PAAD_G12V_HALLMARK_HEDGEHOG_SIGNALING"
    ) +
    theme_fig18_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        title = "Hedgehog signaling",
        tag = "d"
    )

panel_D <- read_fig_proto(
        "rankplot_PAAD_G12V_HALLMARK_HEDGEHOG_SIGNALING"
    ) +
    theme_fig18() +
    theme(
        plot.title = element_text(size = 6, family = "Arial"),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        panel.grid = element_blank()
    )

panel_BCD_legend <- wrap_elements(plot = cowplot::get_legend(panel_D)) +
    theme_fig18() +
    theme(
        plot.margin = margin(0, 0, 0, 0, "mm")
    )

panel_D <- panel_D +
    theme_fig18() +
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


#### ---- E. Genetic dependency heatmap ---- ####
# A heatmap of the genes found to have allele-specific genetic dependencies.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_E <- read_fig_proto("PAAD_CRISPR_manhattan_ward.D2_pheatmap")[[4]]

pre_panel_E_main <- gtable::gtable_filter(pre_panel_E,
                                          "legend",
                                          invert = TRUE)

panel_E <- wrap_elements(plot = pre_panel_E_main) *
    theme_fig18(tag_margin = margin(0, -8, 0, 6, "mm")) %+replace%
    theme(
        plot.margin = margin(-10, -6, -1, -4, "mm")
    ) +
    labs(tag = "e")



prep_pheatmap_legend <- function(name) {
    read_fig_proto(name) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_fig18() +
        theme(
            plot.margin = margin(0, 3, 0, 3, "mm"),
            plot.title = element_text(size = 5, family = "Arial"),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 5, family = "Arial"),
            plot.background = element_rect(fill = NA, color = NA),
            panel.background = element_rect(fill = NA, color = NA)
        )
}

prep_pheatmap_colorbar <- function(name) {
    read_fig_proto(name) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_fig18() +
            theme(
                plot.margin = margin(0, 3, 0, 3, "mm"),
                plot.title = element_text(size = 5, family = "Arial"),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 5, family = "Arial"),
                plot.background = element_rect(fill = NA, color = NA),
                panel.background = element_rect(fill = NA, color = NA)
            )
}

panel_E_legend1 <- prep_pheatmap_colorbar(
        "PAAD_CRISPR_manhattan_ward.D2_pheatmap_heatpal"
    ) +
    labs(title = "scaled\ndep. score")

panel_E_legend2 <- prep_pheatmap_legend(
        "PAAD_CRISPR_manhattan_ward.D2_pheatmap_allelepal"
    ) +
    labs(title = "allele")

panel_E_legend3 <- prep_pheatmap_legend(
        "PAAD_CRISPR_manhattan_ward.D2_pheatmap_clusterpal"
    ) +
    labs(title = "cluster")

panel_E_legend <- (
        panel_E_legend1 | panel_E_legend2 | panel_E_legend3
    ) +
        plot_layout(widths = c(1, 1, 1))

panel_E_legend <- wrap_elements(full = panel_E_legend) +
    theme_fig18()


#### ---- F. Genetic dependency box-plots ---- ####
# Box-plots of example genes found to have allele-specific genetic dependencies.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

panel_F_files <- c(
    "PAAD-CEP350",
    "PAAD-EGLN2",
    "PAAD-JUN",
    "PAAD-NUMB",
    "PAAD-TMED2"
)

panel_F_plots <- as.list(rep(NA, length(panel_F_files)))
names(panel_F_plots) <- panel_F_files
for (f in panel_F_files) {

    panel_F_plots[[f]] <- read_fig_proto(f) +
        theme_fig18(tag_margin = margin(0, 0, 0, -2, "mm")) %+replace%
        theme(
            plot.title = element_text(size = 6, face = "bold"),
            axis.title.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(0, 0, 3, 0, "mm")
        ) +
        labs(
            title = str_remove(file_sans_ext(f), "PAAD-")
        )

}

panel_F <- wrap_plots(panel_F_plots, ncol = 1)
panel_F[[1]] <- panel_F[[1]] + labs(tag = "f")


#### ---- Figure assembly ---- ####


col2_layout <- "
AAAAAAAA
AAAAAAAA
AAAAAAAA
AAAAAAAA
AAAAAAAA
AAAAAAAA
AAAAAAAA
AAAAAAAA
##BBBB##
"

{
    set.seed(0)

    panels_BCD <- (panel_B_density / panel_B /
                   panel_C_density / panel_C /
                   panel_D_density / panel_D) +
        plot_layout(heights = c(2, 5, 2, 5, 2, 5))

    column1 <- (panel_A / wrap_elements(full = panels_BCD) / panel_BCD_legend) +
        plot_layout(heights = c(400, 900, 1))

    column2 <- (panel_E + panel_E_legend) + plot_layout(design = col2_layout)

    # COMPLETE FIGURE
    full_figure <- (column1 | column2 | panel_F) +
        plot_layout(widths = c(1, 3, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
