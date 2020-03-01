# Figure 4. COAD allele-specific dependencies.

FIGNUM <- 4
SUPPLEMENTAL <- FALSE
VERSION <- 1

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig4 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. GSEA ---- ####
# GSEA of the dependencies.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_A <- read_fig_proto("gsea-results-COAD-select.rds", FIGNUM) +
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
            title = "-*log*<sub>10]</sub>(adj. p-value)",
            title.position = "left",
            label.position = "right",
            label.hjust = 0,
            keywidth = unit(0, "mm"),
            order = 2
        )
    ) +
    theme_fig4(tag_margin = margin(0, 0, 0, -3, "mm")) +
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

theme_fig4_densityplots <- function(tag_margin_l = -3) {
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
        "rankline_COAD_G12V_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT.rds",
        FIGNUM
    ) +
    theme_fig4_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        tag = "b",
        title = "Respiratory electron transport"
    )

panel_B <- read_fig_proto(
        "rankplot_COAD_G12V_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT.rds",
        FIGNUM
    ) +
    theme_fig4() +
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
        "rankplot_COAD_G13D_REACTOME_COMPLEMENT_CASCADE.rds",
        FIGNUM
    ) +
    theme_fig4() +
    theme(
        plot.title = element_text(size = 6, family = "Arial"),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        panel.grid = element_blank()
    )

panel_BC_legend <- wrap_elements(plot = cowplot::get_legend(panel_C)) +
    theme_fig4() +
    theme(
        plot.margin = margin(0, 0, 0, 0, "mm")
    )

panel_C_density <- read_fig_proto(
        "rankline_COAD_G13D_REACTOME_COMPLEMENT_CASCADE.rds",
        FIGNUM
    ) +
    theme_fig4_densityplots(tag_margin_l = -7) +
    labs(
        y = "density",
        title = "Complement cascade",
        tag = "c"
    )

panel_C <- panel_C +
    theme_fig4() +
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



#### ---- D. Heatmap of linear model ---- ####
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal in COAD.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_D <- read_fig_proto(
        "COAD_CRISPR_manhattan_ward.D2_pheatmap.rds",
        FIGNUM
    )[[4]]

pre_panel_D_main <- gtable::gtable_filter(pre_panel_D,
                                          "legend",
                                          invert = TRUE)

panel_D <- wrap_elements(plot = pre_panel_D_main) *
    theme_fig4(tag_margin = margin(0, -8, 0, 6, "mm")) %+replace%
    theme(
        plot.margin = margin(-10, -9, -13, -7, "mm")
    ) +
    labs(tag = "d")


prep_pheatmap_legend <- function(name) {
    read_fig_proto(name, FIGNUM) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_fig4() +
        theme(
            plot.title = element_text(size = 5, family = "Arial"),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 5, family = "Arial"),
            plot.background = element_rect(fill = NA, color = NA),
            panel.background = element_rect(fill = NA, color = NA)
        )
}

prep_pheatmap_colorbar <- function(name) {
    read_fig_proto(name, FIGNUM) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_fig4() +
            theme(
                plot.title = element_text(size = 5, family = "Arial"),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 5, family = "Arial"),
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
        x_axis_text <- NULL
    }

    panel_E_plots[[f]] <- read_fig_proto(f, FIGNUM) +
        theme_fig4(tag_margin = margin(0, 0, 0, -2, "mm")) %+replace%
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



################################################################################
#----------------------------------  LUAD  ------------------------------------#
################################################################################



#### ---- F, G. Ranked heatmaps of GSEA ---- ####
# Two heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_F_density <- read_fig_proto("rankline_LUAD_G12C_BIOCARTA_ERK_PATHWAY",
                                  FIGNUM, supp = SUPPLEMENTAL) +
    theme_fig4_densityplots() +
    labs(tag = "f", title = "Erk pathway")

panel_F <- read_fig_proto("rankplot_LUAD_G12C_BIOCARTA_ERK_PATHWAY",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_fig4() +
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


panel_G_density <- read_fig_proto("rankline_LUAD_G12C_PID_BARD1_PATHWAY",
                                  FIGNUM, supp = SUPPLEMENTAL) +
    theme_fig4_densityplots() +
    labs(tag = "g", title = "Bard1 pathway")

panel_G <- read_fig_proto("rankplot_LUAD_G12C_PID_BARD1_PATHWAY",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_fig4() +
    theme(
        plot.title = element_text(size = 6, family = "Arial"),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        panel.grid = element_blank()
    )

panel_FG_legend <- wrap_elements(plot = cowplot::get_legend(panel_G)) +
    theme_fig4() +
    theme(
        plot.margin = margin(1, 0, 0, 0, "mm")
    )

panel_G <- panel_G +
    theme_fig4() +
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



#### ---- H. Heatmap of linear model ---- ####
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal in LUAD.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal in COAD.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_H <- read_fig_proto(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap.rds",
        FIGNUM
    )[[4]]

pre_panel_H_main <- gtable::gtable_filter(pre_panel_H,
                                          "legend",
                                          invert = TRUE)

panel_H <- wrap_elements(plot = pre_panel_H_main) *
    theme_fig4(tag_margin = margin(0, -7, 0, 3, "mm")) +
    theme(
        plot.margin = margin(-10, -9, -13, -2, "mm")
    ) +
    labs(tag = "h")

panel_H_legend1 <- prep_pheatmap_colorbar(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap_heatpal.rds"
    ) +
    labs(title = "scaled\ndep. score")

panel_H_legend2 <- prep_pheatmap_legend(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap_allelepal.rds"
    ) +
    labs(title = "allele")

panel_H_legend3 <- prep_pheatmap_legend(
        "LUAD_CRISPR_manhattan_ward.D2_pheatmap_clusterpal.rds"
    ) +
    labs(title = "cluster")


panel_H_legend <- (
        panel_H_legend1 / panel_H_legend2 / panel_H_legend3 / plot_spacer()
    ) +
        plot_layout(heights = c(1, 1, 1, 3))

panel_H_legend <- wrap_elements(plot = panel_H_legend) +
    theme_fig4() %+replace%
    theme(
        plot.margin = margin(0, 0, 0, 0, "mm")
    )



#### ---- I. PPI of genes from cluster in heatmap ---- ####
# A large PPI composed of genes in cluster 4 of the heatmap (panel h).
# original script: "src/10_15_linear-modeling-syn-let_ppi-subnetworks.R"

panel_I <- read_fig_proto("LUAD_cluster-1_component-1",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_graph_comutation()

panel_I <- wrap_elements(plot = panel_I) +
    theme_graph_comutation() %+replace%
    theme(plot.margin = margin(3, 0, -1, 0, "mm")) +
    labs(tag = "i")


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    ## COAD ##

    panel_BC <- (panel_B_density / panel_B / panel_C_density / panel_C) +
        plot_layout(heights = c(2, 5, 2, 5))

    column_3 <- (
        panel_E /
        wrap_elements(full = panel_D_legend)
    )

    # COMPLETE FIGURE
    coad_figure <- (
        (
            panel_A / wrap_elements(full = panel_BC) / panel_BC_legend +
                plot_layout(heights = c(400, 600, 1))
        ) |
        (
            panel_D
        ) |
        (
            column_3
        )
    ) +
        plot_layout(widths = c(3, 7, 2.1))


    ## LUAD ##

    panel_FG <- (panel_F_density / panel_F / panel_G_density / panel_G)

    luad_figure <- (
        (
            (panel_FG / panel_FG_legend / panel_I) +
            plot_layout(heights = c(200, 500, 200, 500, 1, 800))
        ) |
        (
            panel_H
        ) |
        (
            panel_H_legend
        )
    ) +
        plot_layout(widths = c(10, 12, 2.1))


    ## Full Figure ##
    full_figure <- wrap_elements(full = coad_figure) /
        wrap_elements(full = luad_figure) +
        plot_layout(heights = c(1, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
