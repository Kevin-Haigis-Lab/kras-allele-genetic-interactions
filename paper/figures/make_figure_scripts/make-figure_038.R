# Figure 038. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 38

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 110


theme_fig38 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


# Special theme for graphs from 'ggraph'.
theme_graph_fig38 <- function() {
  theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
      plot.title = element_blank(),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(0, 0, 0, 0, "mm")
      ),
      legend.margin = margin(0, 0, 0, 0, "mm"),
      legend.position = "right",
      legend.title = element_text(size = 6, hjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0),
      plot.margin = margin(0, 0, 0, 0, "mm")
    )
}


#### ---- A. PPIN for PAAD ---- ####
# Annotated PPIN for PAAD.
# original script: "src/40_17_comparing-PAAD-allele-subnetworks.R"

panel_A <- read_fig_proto("paad_overlap_comparison_plot") +
  theme_graph_fig38() +
  theme(
    legend.position = c(0.1, 0.80),
    legend.title = element_blank()
  )


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
