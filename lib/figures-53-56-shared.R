# Functions shared by Figures 53-56

cancer_labels <- mutsig_number_of_samples_labels %>%
  deframe()

make_cancer_label <- function(cancer, tag) {
  tibble(x = 0, y = 0, label = cancer_labels[[cancer]]) %>%
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
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(-1, -1, -1, -1, "mm")
      )
    ) +
    labs(tag = tag)
}
