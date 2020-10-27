# Figure 009. Distribution of mutation burden per sample and types
# of mutations.

FIGNUM <- 9

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig9 <- function() {
  theme_comutation() %+replace%
    theme(
      legend.title = element_blank(),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(-1, -1, -1, -1, "mm")
      )
    )
}


cancer_panel_letters <- letters[1:4]
names(cancer_panel_letters) <- c("COAD", "LUAD", "MM", "PAAD")


# Get the separate pieces for a cancer for the panels in this figure.
get_panel_pieces <- function(cancer) {
  a <- read_fig_proto(
    as.character(glue("{cancer}_coding_muts_distribution"))
  ) +
    theme_fig9() +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_markdown(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text = element_blank()
    ) +
    labs(
      tag = cancer_panel_letters[cancer]
    )
  b <- read_fig_proto(as.character(glue("{cancer}_mutation_types"))) +
    theme_fig9() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if (cancer != "PAAD") {
    b <- b +
      theme(legend.position = "none")
  }

  return(list(a, b))
}


#### ---- A-D. Distribution of mutation burden per sample ---- ####

# The distribution of mutation burden per sample.
# original script: "src/90_03_mutation-burden-distribution.R"

{
  panel_A <- get_panel_pieces("COAD") %>%
    wrap_plots(ncol = 1)

  panel_B <- get_panel_pieces("LUAD") %>%
    wrap_plots(ncol = 1)

  panel_C <- get_panel_pieces("MM") %>%
    wrap_plots(ncol = 1)

  panel_D <- get_panel_pieces("PAAD") %>%
    wrap_plots(ncol = 1)
}



#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (
    (panel_A | panel_B) / (panel_C | panel_D) / guide_area()
  ) +
    plot_layout(
      heights = c(8, 8, 1),
      guides = "collect"
    )

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
