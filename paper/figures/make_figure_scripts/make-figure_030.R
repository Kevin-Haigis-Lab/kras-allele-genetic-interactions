# Figure 030. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 30

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")
FIG_DIMENSIONS$height <- 50


theme_fig30 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}

theme_fig30_densityplots <- function(tag_margin_l = -3) {
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



#### ---- A & B. Ranked heatmaps of GSEA ---- ####
# Two heatmaps showing the ranks of genes in the enriched genesets.
# original script: "src/10_37_gsea-depmap-analysis.R"

x_label <- expression("" %<-% "greater dep. - ranked gene effect - less dep." %->% "")


panel_A_density <- read_fig_proto(
        "rankline_LUAD_G12C_BIOCARTA_P53HYPOXIA_PATHWAY"
    ) +
    theme_fig30_densityplots() +
    labs(tag = "a", title = "p53 hypoxia pathway")

panel_A <- read_fig_proto("rankplot_LUAD_G12C_BIOCARTA_P53HYPOXIA_PATHWAY") +
    theme_fig30() +
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


panel_B_density <- read_fig_proto("rankline_LUAD_G12C_PID_BARD1_PATHWAY") +
    theme_fig30_densityplots() +
    labs(tag = "b", title = "Bard1 pathway")

panel_B <- read_fig_proto("rankplot_LUAD_G12C_PID_BARD1_PATHWAY") +
    theme_fig30() +
    theme(
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "none",
        panel.grid = element_blank()
    )

panel_B <- panel_B +
    theme_fig30() +
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


panel_AB_legend <- tibble(lbl = c("G12C", "G12V", "WT")) %>%
    mutate(lbl = fct_rev(factor_alleles(lbl))) %>%
    ggplot(aes(x = 1, y = lbl, label = lbl, color = lbl)) +
    geom_text(size = 2, fontface = "bold", family = "Arial") +
    scale_color_manual(values = short_allele_pal) +
    theme_void() +
    theme(
        legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold", size = 6)
    ) +
    labs(title = "*KRAS*<br>allele")


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (
        (panel_A_density / panel_A) + plot_layout(heights = c(200, 500)) |
        (panel_B_density / panel_B) + plot_layout(heights = c(200, 500)) |
        (plot_spacer() / panel_AB_legend / plot_spacer())
    ) +
        plot_layout(widths = c(10, 10, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
