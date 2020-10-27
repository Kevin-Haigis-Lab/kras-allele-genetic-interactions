# Figure 050. #> VAF of KRAS and comutation mutations

FIGNUM <- 50

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 120

theme_fig50 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


#### ---- A. KRAS mutation VAF ---- ####
# Distribution of VAF of KRAS mutations.
# original script: "src/90_03_mutation-burden-distribution.R"

panel_A <- read_fig_proto("kras-adj-vaf-distribution") +
  theme_fig50() +
  theme(
    panel.spacing.x = unit(4, "mm"),
    axis.title.x = element_markdown(),
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold")
  ) +
  labs(tag = "a")



#### ---- B. KRAS mutation VAF ---- ####
# Distribution of VAF of KRAS mutations.
# original script: "src/90_03_mutation-burden-distribution.R"

panel_B <- read_fig_proto("comutation-genes-vaf-dist") +
  theme_fig50() +
  theme(
    panel.spacing.x = unit(4, "mm"),
    axis.title.x = element_markdown(),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  labs(tag = "b")


#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (panel_A / panel_B) +
    plot_layout()

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
