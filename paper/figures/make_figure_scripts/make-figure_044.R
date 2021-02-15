# Figure 044. Observed vs. predicted KRAS allele frequencies with all alleles.

FIGNUM <- 44

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 165


theme_fig44 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


#### ---- A. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/50_12_observed-predicted-kras-alleles_v3.R"

panel_A_plots <- paste0(
  c("COAD", "LUAD", "MM", "PAAD"),
  "_predict-ALL-allele-freq_scatter"
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

get_max_val <- function(gg) {
  data <- ggplot_build(read_fig_proto(gg))$plot$data
  return(max(c(data$upper_ci, data$observed_allele_frequency)))
}

panel_A_proto_list <- imap(panel_A_plots, function(x, idx) {
  max_val <- get_max_val(x)
  min_val <- (-0.1 / 0.45) * max_val

  if (idx == 1) {
    color_guide <- panel_A_guide_legend(10)
    shape_guide <- panel_A_guide_legend(20)
  } else {
    color_guide <- FALSE
    shape_guide <- FALSE
  }

  p <- read_fig_proto(x) +
    scale_x_continuous(
      limits = c(min_val, max_val),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_y_continuous(
      limits = c(min_val, max_val),
      expand = expansion(mult = c(0, 0.03)),
      labels = function(x) {
        ifelse(x == 0, "", x)
      }
    ) +
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

  if (idx %in% c(1, 3)) {
    p <- p + labs(y = "observed frequency")
  }

  if (idx %in% c(3, 4)) {
    p <- p + labs(x = "predicted frequency")
  }

  return(p)
})

panel_A <- wrap_plots(panel_A_proto_list, nrow = 2, guides = "collect") &
  theme_fig44() %+replace%
    theme(
      plot.title = element_text(size = 7, vjust = 2),
      plot.subtitle = element_markdown(hjust = 0, vjust = 0.3),
      legend.position = "right",
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5),
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.spacing.y = unit(3, "mm"),
      plot.margin = margin(3, 3, 3, 3, "mm")
    )

panel_A <- (panel_A | guide_area()) +
  plot_layout(widths = c(20, 1))


#### ---- Figure assembly ---- ####

{
  set.seed(0)

  # COMPLETE FIGURE
  full_figure <- panel_A +
    plot_layout()

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
