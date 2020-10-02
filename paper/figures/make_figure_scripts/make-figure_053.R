# Figure 053. #> Distributions and statistics for levels of mutational
#  signatures in tumor samples of various KRAS alleles.

FIGNUM <- 53

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig53 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_comutation() %+replace%
  theme(
    legend.title = element_blank(),
    plot.tag = element_text(size = 7,
                            face = "bold",
                            margin = tag_margin)
  )
}



prepare_ridge_patch <- function(p) {
  p
}


all_panels <- read_fig_proto("ggridge-stats_all-plots") %>%
  map(prepare_ridge_patch)


#### ---- Figure assembly ---- ####
# Rscript: "src/50_35_mutational-signature-allele-associations.R"

make_cancer_label <- function(cancer) {
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
    theme_void(base_size = 7, base_family = "Arial")
}

{
  # COMPLETE FIGURE
  coad_design <- "
  AABB
  CCDD
  "

  mm_design <- "
  AABBCC
  "

  paad_design <- "
  AABBCC
  DDDEEE
  "

  coad_panels <- wrap_plots(all_panels[1:4], design = coad_design)
  luad_panels <- wrap_plots(all_panels[5:8], design = coad_design)
  mm_panels <- wrap_plots(all_panels[9:11], design = mm_design)
  paad_panels <- wrap_plots(all_panels[12:16], design = paad_design)

  full_figure <- (
    ((make_cancer_label("COAD") | coad_panels) + plot_layout(widths = c(1, 30))) /
    ((make_cancer_label("LUAD") | luad_panels) + plot_layout(widths = c(1, 30))) /
    ((make_cancer_label("MM") | mm_panels) + plot_layout(widths = c(1, 30))) /
    ((make_cancer_label("PAAD") | paad_panels) + plot_layout(widths = c(1, 30)))
  ) +
  plot_layout(heights = c(2, 2, 1, 2))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
