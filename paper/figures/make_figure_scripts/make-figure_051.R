# Figure 051. #> Observed vs. predicted cancer-specific KRAS alleles and
# probability of alleles in tumor samples.

FIGNUM <- 51

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 100

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
      legend.text = element_text(size = 6),
      strip.text = element_text(size = 7),
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

codon_pal <- codon_palette[!names(codon_palette) %in% c("Other")]

panel_A_guide_legend <- function(order) {
  guide_legend(
    order = order,
    title.position = "top",
    title.theme = element_markdown(
      hjust = 0.5,
      face = "plain",
      size = 6,
      family = "Arial"
    ),
    keyheight = unit(3, "mm"),
    label.position = "right"
  )
}


panel_A_proto_list <- imap(panel_A_plots, function(x, idx) {
  if (idx == 1) {
    color_guide <- panel_A_guide_legend(10)
    shape_guide <- panel_A_guide_legend(20)
  } else {
    color_guide <- FALSE
    shape_guide <- FALSE
  }
  p <- read_fig_proto(x) +
    scale_color_manual(
      values = codon_pal,
      drop = FALSE,
      breaks = kras_hotspot_codons$char,
      guide = color_guide
    ) +
    scale_shape_manual(
      values = c(17, 16),
      drop = FALSE,
      guide = shape_guide
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
      plot.title = element_text(size = 7, vjust = -1),
      plot.subtitle = element_text(hjust = 0, vjust = -3),
      plot.margin = margin(0, 1, 0, 1, "mm"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.key.width = unit(2, "mm"),
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


#### ---- B. Probability of alleles ---- ####
# Average probability of the KRAS alleles in all tumor samples.
# original script: "src/50_13_observed-predicted-kras-alleles_v3_per_tsb.R"

panel_B <- read_fig_proto("allele_prob_per_allele_plot") +
  theme_fig51() +
  theme(
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.justification = "left",
    legend.box.margin = margin(0, -3, 0, -2, "mm"),
    legend.text = element_markdown()
  ) +
  labs(tag = "b")

#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (panel_A / panel_B) +
    plot_layout(heights = c(3, 2))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
