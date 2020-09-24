# Figure 028. #> Mutational signature main figure.
# (This is a copy of make-fig 28, but without the panel D on mut sig spectra.)

FIGNUM <- 40

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 150

theme_fig40 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


#### ---- A. KRAS hot-spot mutation frequency barplot ---- ####
# The distribution of mutations at the 4 hotspots.
# original script: "src/90_05_kras-allele-distribution.R"

codons_to_label <- c(12, 13, 61, 146)

panel_A <- read_fig_proto("lollipop-kras_hotspot-only") +
  theme_fig40() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.y = element_markdown(),
    legend.position = c(0.80, 0.82),
    legend.background = element_blank(),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    plot.margin = margin(0, 2, 0, 0, "mm")
  ) +
  labs(
    tag = "a",
    y = "number of samples (*log*<sub>10</sub>)"
  )


#### ---- B. Distribution of KRAS alleles ---- ####
# The distribution of KRAS alleles across cancers and codon.
# original script: "src/90_05_kras-allele-distribution.R"

panel_B <- read_fig_proto("allele_dist_dotplot") +
  theme_fig40() +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_markdown(),
    axis.title.y = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold")
  ) +
  labs(tag = "b")


#### ---- C. Percent of samples with KRAS mutation ---- ####
# Percent of samples with a KRAS mutation per cancer.
# original script: "src/90_05_kras-allele-distribution.R"

panel_C <- read_fig_proto("cancer_freq_kras_mut_column") +
  scale_x_continuous(
    expand = expansion(add = c(0, 0.02)),
    breaks = c(0.2, 0.4, 0.6, 0.8),
    labels = function(x) {
      paste0(round(x * 100), "%")
    }
  ) +
  theme_fig40() +
  theme(
    axis.title.x = element_markdown(),
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  labs(tag = "c")


#### ---- D. Mutational signatures probability of causing KRAS allele ---- ####
# The probability that each allele was caused by each detectable mutational
# signature.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

style_mutsig_prob_barplots <- function(plt, i, tag = NULL, y = NULL) {
  themed_plt <- plt +
    scale_fill_manual(
      values = mutsig_descrpt_pal,
      guide = guide_legend(
        nrow = 1,
        title.position = "left",
        title.vjust = 0.2,
        label.position = "top",
        label.hjust = 0.5,
        label.vjust = -4.5
      )
    ) +
    theme_fig40() +
    theme(
      plot.title = element_text(vjust = 0.5, size = 7, face = "bold"),
      plot.margin = margin(1, 1, 1, 1, "mm"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 5.8),
      legend.position = "none",
      legend.title = element_text(size = 6, face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm"),
      legend.spacing.x = unit(1, "mm"),
      legend.spacing.y = unit(0, "mm")
    )
  if (i == 1) {
    themed_plt <- themed_plt + labs(tag = tag)
  }
  if (i %% 2 == 1) {
    themed_plt <- themed_plt + labs(y = y)
  } else {
    themed_plt <- themed_plt + labs(y = NULL)
  }
  return(themed_plt)
}


panel_D <- read_fig_proto("probability-mutsig-caused-allele_barplot-list") %>%
  imap(style_mutsig_prob_barplots, tag = "d", y = "probability")

panel_D_design <- "
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
"

pull_signatures_from_panel_D <- function(x) {
  unique(ggplot_build(x)$plot$data$description)
}

signatures <- map(panel_D, pull_signatures_from_panel_D) %>%
  unlist() %>%
  unique() %>%
  sort() %>%
  as.character()

panel_D_legend <- custom_label_legend(
  signatures,
  y_value = "signature",
  family = "Arial", size = 1.8,
  label.padding = unit(1, "mm"),
  label.size = unit(0, "mm"),
  hjust = 0.5
) +
  scale_fill_manual(values = mutsig_descrpt_pal) +
  theme(
    legend.position = "none",
    plot.margin = margin(-2, 0, -4, 0, "mm"),
    axis.text.y = element_text(size = 6, face = "bold")
  )

panel_D <- wrap_plots(panel_D) +
  plot_layout(design = panel_D_design)


#### ---- E. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/50_12_observed-predicted-kras-alleles_v3.R"

panel_E_plots <- paste0(
  c("COAD", "LUAD", "MM", "PAAD"),
  "_predict-allele-freq_scatter.svg"
)

panel_E_proto_list <- lapply(panel_E_plots, read_fig_proto)
panel_E <- wrap_plots(panel_E_proto_list, nrow = 2, guides = "collect") &
  theme_fig40(tag_margin = margin(0, 0, -3, 0, "mm")) %+replace%
    theme(
      plot.title = element_text(size = 7, face = "bold"),
      plot.margin = margin(1, 1, 1, 1, "mm"),
      legend.position = "bottom",
      legend.title = element_markdown(
        vjust = 0.5, hjust = 0.5,
        face = "bold", size = 6,
        family = "Arial"
      ),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5),
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.spacing.y = unit(3, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(-1, 0, -9, 0, "mm")
    )

for (i in c(2, 4)) panel_E[[i]] <- panel_E[[i]] + labs(y = NULL)
for (i in c(1, 2)) panel_E[[i]] <- panel_E[[i]] + labs(x = NULL)

for (i in seq(1, 4)) {
  panel_E[[i]] <- panel_E[[i]] +
    scale_shape_manual(
      values = c(17, 16),
      drop = FALSE,
      guide = guide_legend(
        title.position = "left",
        label.position = "left"
      )
    )
}

panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")


panel_E <- (panel_E / guide_area()) +
  plot_layout(heights = c(200, 1))


#### ---- Figure assembly ---- ####

{
  set.seed(0)

  # COMPLETE FIGURE
  row_1 <- (panel_A | panel_B | panel_C) +
    plot_layout(widths = c(4, 10, 3))

  panel_D_legend_spacer <- (plot_spacer() | panel_D_legend | plot_spacer()) +
    plot_layout(widths = c(1, 20, 1))

  panel_D_group <- (panel_D / wrap_elements(full = panel_D_legend)) +
    plot_layout(heights = c(15, 1))

  row_2 <- (
    wrap_elements(full = panel_D_group) |
      wrap_elements(full = panel_E)
  ) +
    plot_layout(widths = c(20, 21))

  full_figure <- (row_1 / row_2) +
    plot_layout(heights = c(4, 7))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
