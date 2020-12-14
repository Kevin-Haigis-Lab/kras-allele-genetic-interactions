# Figure 033. GSEA of DepMap of PAAD cell lines

FIGNUM <- 33

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(1, "tall")


theme_fig33 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
  theme_comutation() %+replace%
    theme(
      legend.title = element_blank(),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      )
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
      title = "-log<sub>10</sub>(adj. p-value)",
      title.position = "left",
      label.position = "right",
      label.hjust = 0,
      keywidth = unit(0, "mm"),
      order = 2
    )
  ) +
  theme_fig33(tag_margin = margin(0, 0, 0, -3, "mm")) +
  theme(
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(-2, 0, 0, 0, "mm"),
    legend.title = element_markdown(
      size = 6, angle = 90,
      vjust = 0.5, hjust = 0.5
    ),
    legend.spacing = unit(0, "mm")
  ) +
  labs(tag = "a")


#### ---- B. GSEA ranked-heatmap (1) ---- ####
# A heatmap of an enriched gene set with the cell lines ranked by their
# dependency score for each gene. A density plot along the top helps
# highlight the trend.
# original script: "src/10_37_gsea-depmap-analysis.R"


theme_fig33_densityplots <- function(tag_margin_l = -3) {
  theme_classic_comutation() +
    theme(
      plot.tag = element_text(
        size = 7,
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

x_label <- "← greater dep. - ranked by dep. score - less dep. →"

panel_B_density <- read_fig_proto(
  "rankline_PAAD_G12D_REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT"
) +
  theme_fig33_densityplots(tag_margin_l = -7) +
  labs(
    y = "density",
    tag = "b",
    title = "G2/M DNA damage checkpoint"
  )

panel_B <- read_fig_proto(
  "rankplot_PAAD_G12D_REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT"
) +
  theme_fig33() +
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
  "rankline_PAAD_G12R_REACTOME_PI_3K_CASCADE:FGFR1"
) +
  theme_fig33_densityplots(tag_margin_l = -7) +
  labs(
    y = "density",
    tag = "c",
    title = "PI3K-FGFR1 cascade"
  )

panel_C <- read_fig_proto(
  "rankplot_PAAD_G12R_REACTOME_PI_3K_CASCADE:FGFR1"
) +
  theme_fig33() +
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
  "rankline_PAAD_G12V_REACTOME_CELLULAR_SENESCENCE"
) +
  theme_fig33_densityplots(tag_margin_l = -7) +
  labs(
    y = "density",
    title = "Cellular senescence",
    tag = "d"
  )

panel_D <- read_fig_proto(
  "rankplot_PAAD_G12V_REACTOME_CELLULAR_SENESCENCE"
) +
  theme_fig33() +
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


lbl_alleles <- ggplot_build(panel_D)$plot$data$kras_allele %>%
  unique() %>%
  factor_alleles() %>%
  sort() %>%
  as.character()

panel_BCD_legend <- custom_label_legend(
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
    axis.text.y = element_markdown(
      hjust = 0.5, face = "bold",
      size = 6, family = "Arial"
    )
  )



#### ---- Figure assembly ---- ####

{
  set.seed(0)

  panels_BCD <- (panel_B_density / panel_B /
    panel_C_density / panel_C /
    panel_D_density / panel_D) +
    plot_layout(heights = c(2, 5, 2, 5, 2, 5))

  # COMPLETE FIGURE
  full_figure <- (panel_A / wrap_elements(full = panels_BCD) / panel_BCD_legend) +
    plot_layout(heights = c(20, 40, 1))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
