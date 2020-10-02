# Figure 055. #> Distributions and statistics for probability of
# mutational signatures to cause KRAS mutations (box-plots).

FIGNUM <- 56

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig56 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_comutation() %+replace%
  theme(
    legend.title = element_blank(),
    plot.tag = element_text(size = 7,
                            face = "bold",
                            margin = tag_margin)
  )
}


prepare_ridge_patch <- function(p, idx) {

  i_to_turn_xaxis_text <- c(1:3, 10:13)
  if (idx %in% i_to_turn_xaxis_text) {
    p <- p +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      panel.grid.major.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  }
  return(p)
}


all_panels <- read_fig_proto("box_plots_sig_cause") %>%
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
  ABC
  "

  luad_design <- "
  ABCD
  "

  mm_design <- "
  AB
  "

  paad_design <- "
  ABCD
  "

  coad_panels <- wrap_plots(all_panels[1:3], design = coad_design)
  luad_panels <- wrap_plots(all_panels[4:7], design = luad_design)
  mm_panels <- wrap_plots(all_panels[8:9], design = mm_design)
  paad_panels <- wrap_plots(all_panels[10:13], design = paad_design)

  full_figure <- (
    (
      (make_cancer_label("COAD", "a") | coad_panels) +
      plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("LUAD", "b") | luad_panels) +
      plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("MM", "c") | mm_panels) +
      plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("PAAD", "d") | paad_panels) +
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
