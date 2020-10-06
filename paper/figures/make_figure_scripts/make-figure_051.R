# Figure 051. #> Observed vs. predicted cancer-specific KRAS alleles and
# probability of alleles in tumor samples.

FIGNUM <- 51

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 140

fig51_tag_element <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  element_text(
    size = 7,
    family = "arial",
    face = "bold",
    margin = tag_margin
  )
}

theme_fig51 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_comutation() %+replace%
    theme(
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 6),
      strip.text = element_text(size = 7, face = "bold"),
      plot.tag = fig51_tag_element(tag_margin)
    )
}


#### ---- A. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/50_12_observed-predicted-kras-alleles_v3.R"

panel_A_plots <- paste0(
  c("COAD", "LUAD", "MM", "PAAD"),
  "_predict-allele-freq_scatter.svg"
)

codon_pal <- codon_palette[names(codon_palette) != "Other"]

panel_A_guide_legend <- function(order) {
  guide_legend(
    order = order,
    title.position = "top",
    title.theme = element_markdown(
      hjust = 0.5,
      face = "bold",
      size = 6,
      family = "Arial"
    ),
    keyheight = unit(3, "mm"),
    label.position = "right"
  )
}


panel_A_proto_list <- map(panel_A_plots, function(x) {
  p <- read_fig_proto(x) +
    scale_color_manual(
      values = codon_pal,
      drop = FALSE,
      guide = panel_A_guide_legend(10)
    ) +
    scale_shape_manual(
      values = c(17, 16),
      drop = FALSE,
      guide = panel_A_guide_legend(20)
    ) +
    labs(
      x = "predicted frequency",
      y = "observed frequency"
    )
  return(p)
})

panel_A <- wrap_plots(panel_A_proto_list, nrow = 1, guides = "collect") &
  theme_fig51() %+replace%
    theme(
      plot.title = element_text(size = 7, vjust = 1, face = "bold"),
      plot.subtitle = element_markdown(hjust = 0, vjust = 0),
      plot.margin = margin(0, 1, 0, 1, "mm"),
      legend.position = "right",
      legend.background = element_blank()
    )

for (i in seq(2, 4)) {
  panel_A[[i]] <- panel_A[[i]] + labs(y = NULL)
}

panel_A <- (panel_A | guide_area()) +
  plot_layout(widths = c(3, 3, 3, 3, 1))
panel_A <- wrap_elements(full = panel_A) +
  labs(tag = "a") +
  theme(
    plot.tag = fig51_tag_element()
  )



#### ---- B. Accuracy to predict allele ---- ####
# Accuracy of the mutational signatures to predict observed allele.
# original script: "src/50_13_observed-predicted-kras-alleles_v3_per_tsb.R"

panel_B <- read_fig_proto("allele_accuracy_barplots") +
  theme_fig51() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.spacing = unit(2, "mm"),
    legend.spacing.y = unit(2, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.title = element_text(face = "bold"),
    legend.justification = "left"
  ) +
  labs(tag = "b")



#### ---- C. Probability of alleles ---- ####
# Average probability of the KRAS alleles in all tumor samples.
# original script: "src/50_13_observed-predicted-kras-alleles_v3_per_tsb.R"

panel_C <- read_fig_proto("allele_prob_per_allele_plot") +
  theme_fig51() +
  theme(
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height = unit(3, "mm"),
    legend.justification = "left"
  ) +
  labs(tag = "c")



#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (
    panel_A /
      panel_B /
      panel_C
  ) +
    plot_layout(heights = c(6, 4, 4))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
