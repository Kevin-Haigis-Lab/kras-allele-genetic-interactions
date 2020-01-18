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
        legend.title = element_blank(),
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
# Heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

x_label <- expression("" %<-% "greater dep. - ranked gene effect - less dep." %->% "")

panel_B <- read_fig_proto(
        "rankplot_COAD_G12V_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT.rds",
        FIGNUM
    ) +
    theme_fig4() +
    theme(
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    labs(
         title = "Respiratory electron transport",
         x = x_label)


#### ---- C. Ranked heatmaps of GSEA ---- ####
# Heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

panel_C <- read_fig_proto(
        "rankplot_COAD_G13D_KEGG_OXIDATIVE_PHOSPHORYLATION.rds",
        FIGNUM
    ) +
    theme_fig4() +
    theme(
        axis.title.y = element_blank(),
        plot.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.direction = "horizontal",
        panel.grid = element_blank()
    ) +
    labs(x = x_label)


#### ---- D. Heatmap of linear model ---- ####
# Heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

pre_panel_D <- read_fig_proto(
        "COAD_CRISPR_manhattan_ward.D2_pheatmap.rds",
        FIGNUM
    )[[4]]

print(pre_panel_D)
pre_panel_D_main <- gtable::gtable_filter(pre_panel_D,
                                          "legend",
                                          invert = TRUE)
pre_panel_D_main$layout
# pre_panel_D_main$layout$b[[10]] <- 2.5
# pre_panel_D_main$heights[4] <- pre_panel_D$heights[1]

panel_D <- wrap_elements(plot = pre_panel_D_main) *
    theme_fig4(tag_margin = margin(-5, -8, 0, 6, "mm")) +
    theme(
        plot.margin = margin(-10, -11, -13, -7, "mm")
    ) +
    labs(tag = "d")


panel_D_legend1_label <- grid::textGrob(
    "scaled dep. score",
    rot = 90,
    gp = grid::gpar(fontsize = 5,
                    fontfamily = "Arial",
                    fontface = "bold"))
panel_D_legend1_label <- wrap_elements(panel = panel_D_legend1_label)

panel_D_legend1 <- pre_panel_D %>%
    gtable::gtable_filter("legend") %>%
    gtable::gtable_filter("annotation", invert = TRUE)
# panel_D_legend1$layout$clip <- "on"
panel_D_legend1 <- wrap_elements(full = panel_D_legend1)

panel_D_legend2 <- pre_panel_D %>%
    gtable::gtable_filter("annotation_legend")
panel_D_legend2$layout$clip <- "on"
panel_D_legend2 <- wrap_elements(full = panel_D_legend2)

panel_D_legend <- panel_D_legend1_label | panel_D_legend1 |
    plot_spacer() |
    panel_D_legend2 &
    theme_fig4()


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    panels_BC <- panel_B / panel_C / guide_area() +
        plot_layout(heights = c(50, 50, 1), guides = "collect")
    panels_BC <- wrap_elements(full = panels_BC) +
        labs(tag = "b") +
        theme_fig4()

    # COMPLETE FIGURE
    full_figure <- (
        (
            panel_A / panels_BC +
            plot_layout(heights = c(2, 3))
        ) |
        (
            panel_D
        ) |
        (
            wrap_elements(full = panel_D_legend) / plot_spacer() +
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