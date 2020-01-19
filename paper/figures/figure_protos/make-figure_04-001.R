# Figure 4. COAD allele-specific dependencies.

FIGNUM <- 4
SUPPLEMENTAL <- FALSE
VERSION <- 1

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")
FIG_DIMENSIONS$height <- FIG_DIMENSIONS$height / 2


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
            barheight = unit(20, "mm"),
            barwidth = unit(2, "mm"),
            order = 1
        )
    ) +
    scale_size_continuous(
        range = c(0.5, 3),
        guide = guide_legend(
            title = expression(paste(italic("log")[10], "(adj. p-val.)")),
            title.position = "left",
            label.position = "right",
            label.hjust = 0,
            keywidth = unit(1, "mm"),
            order = 2
        )
    ) +
    theme_fig4(tag_margin = margin(-5, 0, 0, 0, "mm")) +
    theme(
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        legend.title = element_text(angle = 90),
        legend.spacing = unit(0, "mm")
    ) +
    labs(tag = "a")


#### ---- B. Ranked heatmaps of GSEA ---- ####
# Two heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

x_label <- expression("" %<-% "greater dep. - ranked gene effect - less dep." %->% "")

panel_B1 <- read_fig_proto(
        "rankplot_COAD_G12V_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT.rds",
        FIGNUM
    ) +
    theme_fig4() +
    theme(
        plot.title = element_text(size = 6, family = "Arial"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    labs(title = "Respiratory electron transport",
         x = x_label)

panel_B2 <- read_fig_proto(
        "rankplot_COAD_G13D_REACTOME_COMPLEMENT_CASCADE.rds",
        FIGNUM
    ) +
    theme_fig4() +
    theme(
        plot.title = element_text(size = 7, family = "Arial"),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        panel.grid = element_blank()
    ) +
    labs(title = "Complement cascade",
         x = x_label)


#### ---- c. Heatmap of linear model ---- ####
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal.
# original script: "src/10_10_linear-modeling-syn-let.R"

pre_panel_C <- read_fig_proto(
        "COAD_CRISPR_manhattan_ward.D2_pheatmap.rds",
        FIGNUM
    )[[4]]

pre_panel_C_main <- gtable::gtable_filter(pre_panel_C,
                                          "legend",
                                          invert = TRUE)

panel_C <- wrap_elements(plot = pre_panel_C_main) *
    theme_fig4(tag_margin = margin(-5, -8, 0, 6, "mm")) +
    theme(
        plot.margin = margin(-10, -11, -13, -7, "mm")
    ) +
    labs(tag = "c")


panel_C_legend1_label <- grid::textGrob(
    "scaled dep. score",
    rot = 90,
    gp = grid::gpar(fontsize = 5,
                    fontfamily = "Arial",
                    fontface = "bold"))
panel_C_legend1_label <- wrap_elements(panel = panel_C_legend1_label)

panel_C_legend1 <- pre_panel_C %>%
    gtable::gtable_filter("legend") %>%
    gtable::gtable_filter("annotation", invert = TRUE)
# panel_C_legend1$layout$clip <- "on"
panel_C_legend1 <- wrap_elements(full = panel_C_legend1)

panel_C_legend2 <- pre_panel_C %>%
    gtable::gtable_filter("annotation_legend")
panel_C_legend2$layout$clip <- "on"
panel_C_legend2 <- wrap_elements(full = panel_C_legend2)

panel_C_legend <- panel_C_legend1_label | panel_C_legend1 |
    plot_spacer() |
    panel_C_legend2 &
    theme_fig4()


#### ---- E. Heatmap of linear model ---- ####
# Clustered (pretty) heatmap of genes found to be differentially synthetic
# lethal.
# original script: "src/10_10_linear-modeling-syn-let.R"

panel_D_files <- c(
    "COAD-IDH1.rds",
    "COAD-KNTC1.rds",
    "COAD-PIP5K1A.rds",
    "COAD-WDR26.rds"
)

panel_D_plots <- as.list(rep(NA, length(panel_D_files)))
names(panel_D_plots) <- panel_D_files
for (f in panel_D_files) {
    panel_D_plots[[f]] <- read_fig_proto(f, FIGNUM) +
        theme_fig4() +
        theme(
            axis.title.x = element_blank(),
            legend.position = "none"
        )
}

panel_D <- wrap_plots(panel_D_plots, ncol = 1)
panel_D[[1]] <- panel_D[[1]] + labs(tag = "d")


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    panel_B <- panel_B1 / panel_B2 / guide_area() +
        plot_layout(heights = c(50, 50, 1), guides = "collect")
    panel_B <- wrap_elements(full = panel_B) +
        labs(tag = "b") +
        theme_fig4()

    # COMPLETE FIGURE
    full_figure <- (
        (
            panel_A / panel_B +
            plot_layout(heights = c(2, 3))
        ) |
        (
            panel_C
        ) |
        (
            wrap_elements(full = panel_C_legend) / panel_D +
                plot_layout(heights = c(2, 5))
        )
    ) +
        plot_layout(widths = c(3, 7, 2.1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}


# TODO:
# I am re-running the GSEA of the DepMap data.
# Once it is done, I need to rerun "src/10_37_..."