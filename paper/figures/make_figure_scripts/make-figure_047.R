# Figure 047. #> The first group of results from the analysis for explaining
# allele-specific genetic dependency by comutation (supplementary examples).

FIGNUM <- 47

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 120


theme_fig47 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


theme_coefplot_fig47 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
  theme_minimal_comutation() %+replace%
    theme(
      plot.title = element_markdown(
        size = 7,
        face = "bold",
        family = "Arial",
        hjust = 0
      ),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      ),
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_markdown()
    )
}


theme_lineplot_fig47 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
  theme_fig45(tag_margin = tag_margin) %+replace%
    theme(
      legend.position = "none",
      axis.title.y = element_markdown(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.title.x = element_blank()
    )
}


#### ---- A. SMAD4 explaining ABI1 ---- ####
# SMAD4 explaining dep. of G12D on ABI1 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("PAAD_G12D_ABI1_coef-plot") +
  theme_coefplot_fig47() +
  labs(tag = "a")
panel_A_2 <- read_fig_proto("PAAD_G12D_ABI1_line-plot") +
  theme_lineplot_fig47()



#### ---- B. SMAD4 explaining MYBL2 ---- ####
# SMAD4 explaining dep. of G12D on MYBL2 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12D_MYBL2_coef-plot") +
  theme_coefplot_fig47() +
  labs(tag = "b")
panel_B_2 <- read_fig_proto("PAAD_G12D_MYBL2_line-plot") +
  theme_lineplot_fig47()




#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (
    (panel_A_1 / panel_A_2) + plot_layout(height = c(1, 3)) |
      (panel_B_1 / panel_B_2) + plot_layout(height = c(1, 3))
  )

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
