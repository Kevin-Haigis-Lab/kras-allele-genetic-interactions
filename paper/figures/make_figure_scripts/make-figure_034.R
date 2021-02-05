# Figure 034. Supplemental with gene-level allele-specific dep. for PAAD.

FIGNUM <- 34

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(1, "tall")
FIG_DIMENSIONS$width <- 120

theme_fig34 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
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


#### ---- A. Genetic dependency heatmap ---- ####
# A heatmap of the genes found to have allele-specific genetic dependencies.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_A <- read_fig_proto("PAAD_CRISPR_euclidean_ward.D2_pheatmap")[[4]]

pre_panel_A_main <- gtable::gtable_filter(pre_panel_A,
  "legend",
  invert = TRUE
)

panel_A <- wrap_elements(plot = pre_panel_A_main) *
  theme_fig34(tag_margin = margin(0, -10, 0, 6, "mm")) +
  theme(
    plot.margin = margin(-20, -8, -1, -6, "mm")
  ) +
  labs(tag = "a")


panel_A_legend1 <- read_fig_proto(
  "PAAD_CRISPR_manhattan_ward.D2_pheatmap_heatpal"
)

panel_A_legend <- ggplot_build(panel_A_legend1)$plot$data %>%
  mutate(
    x = 1:n(),
    name = as.character(round(name, 1)),
    name = fct_inorder(name)
  ) %>%
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


#### ---- B. Genetic dependency box-plots ---- ####
# Box-plots of example genes found to have allele-specific genetic dependencies.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

panel_B_files <- c(
  "PAAD-KHDRBS1_extra",
  "PAAD-EGLN2_extra",
  "PAAD-JUN_extra",
  "PAAD-MAPK8_extra",
  "PAAD-BRI3BP_extra"
)

panel_B_plots <- as.list(rep(NA, length(panel_B_files)))
names(panel_B_plots) <- panel_B_files
for (f in panel_B_files) {
  panel_B_plots[[f]] <- read_fig_proto(f) +
    theme_fig34(tag_margin = margin(0, 0, 0, -2, "mm")) %+replace%
    theme(
      plot.title = element_text(size = 6, face = "bold.italic"),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 3, 0, "mm")
    ) +
    labs(
      title = str_remove_all(file_sans_ext(f), "PAAD-|_extra"),
      y = "dependency score"
    )
}

panel_B <- wrap_plots(panel_B_plots, ncol = 1)
panel_B[[1]] <- panel_B[[1]] + labs(tag = "b")



#### ---- Figure assembly ---- ####


{
  set.seed(0)

  # COMPLETE FIGURE
  full_figure <- (
    panel_A |
      (
        (panel_B / wrap_elements(full = panel_A_legend)) +
          plot_layout(heights = c(5, 5, 5, 5, 5, 1))
      )
  ) +
    plot_layout(widths = c(7, 2))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
