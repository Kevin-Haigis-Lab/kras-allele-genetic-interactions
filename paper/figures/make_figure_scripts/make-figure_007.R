# Figure 007. Distribution of mutational signatures across alleles

FIGNUM <- 7

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig7 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
  theme_comutation() %+replace%
    theme(
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      )
    )
}


#### ---- A. Distirubiton of mutational signatures in each sample ---- ####
# The distribution of mutational signature levels in each sample (barplot).
# original script: "src/50_20_mutsignatures-distributions.R"

panel_A <- read_fig_proto("signature-level-per-sample") +
  scale_fill_manual(
    values = mutsig_descrpt_pal,
    guide = guide_legend(
      nrow = 1,
      label.position = "top",
      label.hjust = 0.5,
      label.vjust = -11
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    y = "signature level",
    fill = "signature",
    tag = "a"
  ) +
  theme_fig7() +
  theme(
    legend.position = "none",
    legend.title = element_text(vjust = 0.25),
    axis.text.x = element_blank(),
    legend.key.size = unit(4, "mm"),
    strip.text = element_text(size = 7)
  )


signatures <- ggplot_build(panel_A)$plot$data %>%
  filter(contribution > 0) %>%
  mutate(description = factor(description,
    levels = names(mutsig_descrpt_pal)
  )) %>%
  u_pull(description) %>%
  unlist() %>%
  unique() %>%
  sort() %>%
  as.character()

panel_AB_legend <- custom_label_legend(
  signatures,
  y_value = "signature",
  family = "Arial",
  size = 2.0,
  label.padding = unit(1, "mm"),
  label.size = unit(0, "mm"),
  hjust = 0.5
) +
  scale_fill_manual(values = mutsig_descrpt_pal) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, -5, 0, "mm"),
    axis.text.y = element_text(size = 6, face = "bold")
  )


#### ---- B. Distirubiton of mutational signature levels ---- ####
# The distribution of mutational signature levels in each sample (boxplot).
# original script: "src/50_20_mutsignatures-distributions.R"

panel_B <- read_fig_proto("signature-level-boxplots_with0") +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    tag = "b",
    y = "signature level"
  ) +
  theme_fig7() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 7)
  )


#### ---- D. Levels of clock vs. non-clock mutational signatures ---- ####
# The levels of clock and non-clock mutational signatures per cancer.
# original script: "src/50_20_mutsignatures-distributions.R"

panel_C <- read_fig_proto("clock-signatures_violin-box") +
  facet_wrap(~cancer, nrow = 1) +
  labs(tag = "c") +
  theme_fig7() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.text = element_text(size = 7),
    panel.spacing.x = unit(4, "mm")
  )


#### ---- Figure assembly ---- ####

{
  set.seed(0)

  # COMPLETE FIGURE
  full_figure <- (panel_A + panel_B + plot_layout(widths = c(1))) /
    (
      (
        plot_spacer() | panel_AB_legend | plot_spacer()
      ) +
        plot_layout(widths = c(1, 3, 1))
    ) /
    panel_C +
    plot_layout(heights = c(42, 1, 8))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
