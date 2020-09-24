# Figure 048. #> The group of results from the analysis for explaining
# allele-specific genetic dependency by comutation.

FIGNUM <- 48

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 100

theme_fig48 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


theme_fig48 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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


theme_coefplot_fig48 <- function(tag_margin = margin(-1, -1, -2, -1, "mm")) {
  theme_minimal_comutation() %+replace%
    theme(
      plot.title = element_markdown(
        size = 7,
        face = "bold",
        family = "Arial",
        hjust = 0,
        halign = 0
      ),
      plot.subtitle = element_markdown(
        size = 6,
        face = "bold",
        family = "Arial",
        hjust = 0,
        halign = 0
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


theme_boxplot_fig48 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_fig48(tag_margin = tag_margin) %+replace%
    theme(
      legend.position = "none",
      axis.title.y = element_markdown(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.title.x = element_blank()
    )
}


theme_fig48_void <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_graph_comutation() %+replace%
    theme(
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      ),
      plot.margin = margin(1, -1, -1, -1, "mm")
    )
}

get_coef_plot_title <- function(cancer, gene) {
  glue("{cancer}<br>Effect of mutation on *{gene}* dep.")
}

get_coef_plot_subtitle <- function(gene) {
  glue("Effect of mutation on *{gene}* dep.")
}


#### ---- A. TP53 explaining STARD9 ---- ####
# TP53 explaining dep. of G12D on STARD9 in COAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("COAD_G12D_STARD9_coef-plot") +
  scale_x_continuous(expand = expansion(mult = c(0.17, 0.15))) +
  theme_coefplot_fig48() +
  labs(
    tag = "a",
    x = "estimated effect",
    title = "COAD",
    subtitle = get_coef_plot_subtitle("STARD9")
  )
panel_A_2 <- read_fig_proto("COAD_G12D_STARD9_line-plot") +
  theme_boxplot_fig48()
panel_A_3 <- read_fig_proto("COAD_G12D_TP53_comut-graph-plot.rds") +
  scale_x_continuous(expand = expansion(mult = c(0.50, 0.50))) +
  theme_fig48_void()

panel_A <- (panel_A_1 / panel_A_2) +
  plot_layout(height = c(3, 7))



#### ---- B. SMAD4 explaining EEF1E1 ---- ####
# SMAD4 explaining dep. of G12D on EEF1E1 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12D_EEF1E1_coef-plot") +
  scale_x_continuous(expand = expansion(mult = c(0.17, 0.15))) +
  theme_coefplot_fig48() +
  labs(
    tag = "b",
    x = "estimated effect",
    title = "PAAD",
    subtitle = get_coef_plot_subtitle("EEF1E1")
  )
panel_B_2 <- read_fig_proto("PAAD_G12D_EEF1E1_line-plot") +
  theme_boxplot_fig48()
panel_BCD_3 <- read_fig_proto("PAAD_G12D_SMAD4_comut-graph-plot.rds") +
  scale_x_continuous(expand = expansion(mult = c(0.50, 0.50))) +
  theme_fig48_void()

panel_B <- (panel_B_1 / panel_B_2) +
  plot_layout(height = c(3, 7))


#### ---- C. SMAD4 explaining ABI1 ---- ####
# SMAD4 explaining dep. of G12D on ABI1 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_C_1 <- read_fig_proto("PAAD_G12D_ABI1_coef-plot") +
  scale_x_continuous(expand = expansion(mult = c(0.17, 0.15))) +
  theme_coefplot_fig48() +
  labs(
    tag = "c",
    x = "estimated effect",
    title = "PAAD",
    subtitle = get_coef_plot_subtitle("ABI1")
  )
panel_C_2 <- read_fig_proto("PAAD_G12D_ABI1_line-plot") +
  theme_boxplot_fig48()


panel_C <- (panel_C_1 / panel_C_2) +
  plot_layout(height = c(3, 7))



#### ---- D. SMAD4 explaining MYBL2 ---- ####
# SMAD4 explaining dep. of G12D on MYBL2 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_D_1 <- read_fig_proto("PAAD_G12D_MYBL2_coef-plot") +
  scale_x_continuous(expand = expansion(mult = c(0.17, 0.15))) +
  theme_coefplot_fig48() +
  labs(
    tag = "d",
    x = "estimated effect",
    title = "PAAD",
    subtitle = get_coef_plot_subtitle("MYBL2")
  )
panel_D_2 <- read_fig_proto("PAAD_G12D_MYBL2_line-plot") +
  theme_boxplot_fig48()

panel_D <- (panel_D_1 / panel_D_2) +
  plot_layout(height = c(3, 7))



#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (
    ((
      panel_A | panel_B | panel_C | panel_D
    ) + plot_layout(widths = c(1, 1, 1, 1))) /
      ((
        panel_A_3 | panel_BCD_3
      ) + plot_layout(widths = c(1, 3)))
  ) +
    plot_layout(heights = c(10, 2))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
