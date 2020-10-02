# Figure 055. #> Distributions and statistics for levels of mutational
# signatures in tumors with different KRAS alleles (box-plots).

FIGNUM <- 55

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig55 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_comutation() %+replace%
  theme(
    legend.title = element_blank(),
    plot.tag = element_text(size = 7,
                            face = "bold",
                            margin = tag_margin)
  )
}


prepare_ridge_patch <- function(p, idx) {
  p
}


all_panels <- read_fig_proto("box_plots_sig_levels") %>%
  imap(prepare_ridge_patch)


#### ---- Figure assembly ---- ####

make_cancer_label <- function(cancer, tag) {
  tibble(x = 0, y = 0, label = cancer) %>%
    ggplot(aes(x, y)) +
    geom_text(
      aes(label = label),
      angle = 90,
      vjust = 0,
      hjust = 0.5,
      family = "Arial",
      size = 2.5,
      fontface = "bold"
    ) +
    geom_line(
      data = tibble(x = 1, y = c(-1, 1)),
      color = cancer_palette[[cancer]],
      size = 1
    ) +
    scale_x_continuous(
      limits = c(-1, 2),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    theme_void(base_size = 7, base_family = "Arial") +
    theme(
      plot.tag = element_text(size = 7,
                              face = "bold",
                              margin = margin(-1, -1, -1, -1, "mm"))
    ) +
    labs(tag = tag)
}

{
  # COMPLETE FIGURE
  coad_design <- "
  AABB
  "

  luad_design <- "
  AABB
  "

  mm_design <- "
  AABBCC
  "

  paad_design <- "
  AABBCC
  "

  coad_panels <- wrap_plots(all_panels[1:2], design = coad_design)
  luad_panels <- wrap_plots(all_panels[3:4], design = luad_design)
  mm_panels <- wrap_plots(all_panels[5:7], design = mm_design)
  paad_panels <- wrap_plots(all_panels[8:10], design = paad_design)

  full_figure <- (
    (
      (make_cancer_label("COAD", tag = "a") | coad_panels) +
      plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("LUAD", tag = "b") | luad_panels) +
      plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("MM", tag = "c") | mm_panels) +
      plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("PAAD", tag = "d") | paad_panels) +
      plot_layout(widths = c(1, 30))
    )
  ) +
  plot_layout(heights = c(2, 2, 2, 2))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}