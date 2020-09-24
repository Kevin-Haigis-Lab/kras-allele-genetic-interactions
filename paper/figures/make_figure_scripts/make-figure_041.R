# Figure 041. Heatmap of comutation for MM.

FIGNUM <- 41

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 100

theme_fig41 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


#### ---- A. Labeled MM comutation graph ---- ####
# The comutation graph for MM with every node labeled.
# original script: "src/60_10_MM-specific-oncogenes.R"

panel_A_1 <- read_fig_proto("mm_comut_heatmap_TRUNCATED") +
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "tomato",
    mid = "grey90",
    midpoint = 0.10,
    guide = guide_colorbar(
      barheight = unit(2, "mm")
    )
  ) +
  theme_fig41() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 1, "mm")
  )
panel_A_2 <- read_fig_proto("allele_freq_barplot_TRUNCATED") +
  scale_y_continuous(
    breaks = c(10, 50, 200, 700),
    expand = expansion(mult = c(0, 0.05)),
    trans = "log10"
  ) +
  theme_fig41(margin(-1.9, 0, 0, 0, "mm")) +
  theme(
    plot.title = element_text(
      hjust = 0, face = "bold",
      size = 7, vjust = 5
    ),
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
  labs(y = "num. tumor<br>samples<br>(*log*<sub>10</sub>)")
panel_A_3 <- read_fig_proto("gene_freq_barplot_TRUNCATED") +
  scale_y_continuous(
    breaks = c(25, 100, 200),
    expand = expansion(mult = c(0, 0.05)),
  ) +
  theme_fig41() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 6),
    plot.margin = margin(0, 0, 0, -0.5, "mm")
  )

panel_A_design <- "
    11111#
    222223
    222223
    222223
    222223
"

panel_A <- panel_A_2 + panel_A_1 + panel_A_3 +
  plot_layout(design = panel_A_design)




#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (plot_spacer() | panel_A | plot_spacer()) +
    plot_layout(widths = c(1, 4, 1))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
