# Figure 039. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 39

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 90

theme_fig39 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
  theme_comutation() %+replace%
    theme(
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      )
    )
}

theme_minimal_fig39 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
  theme_minimal_comutation() %+replace%
    theme(
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      )
    )
}


#### ---- A. Lollipop of STK11 for G12C vs rest ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/20_70_luad-g12c-stk11.R"

panel_A <- read_fig_proto("stk11_lollipop_patch")
panel_A[[1]] <- panel_A[[1]] +
  theme_fig39() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(0, 0, -2, 0, "mm"),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = "none"
  )
panel_A[[2]] <- panel_A[[2]] +
  theme_minimal_fig39() +
  theme(
    legend.title = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = "none",
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(-2, 0, -2, 0, "mm")
  )
panel_A[[3]] <- panel_A[[3]] +
  theme_fig39() +
  theme(
    plot.margin = margin(-2, 0, 0, 0, "mm"),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = "none"
  )


variants <- ggplot_build(panel_A[[1]])$plot$data$mutation_type %>%
  unique() %>%
  sort()

panel_A_leg_df <- custom_label_legend_df(variants, colors = "white", )

panel_A_leg <- custom_label_legend_plot(panel_A_leg_df,
  y_value = "mutation type",
  family = "Arial", size = 1.8,
  label.padding = unit(1, "mm"),
  label.size = unit(0, "mm"),
  hjust = 0.5
) +
  scale_fill_manual(values = mod_variant_pal, guide = FALSE) +
  theme(
    axis.text.y = element_text(size = 6, hjust = 0, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_blank()
  )


#### ---- Figure assembly ---- ####

{
  # COMPLETE FIGURE
  full_figure <- (
    (
      (plot_spacer() | panel_A | plot_spacer()) +
        plot_layout(widths = c(1, 3, 1))
    ) /
      (
        (plot_spacer() | panel_A_leg | plot_spacer()) +
          plot_layout(widths = c(1, 2, 1))
      )
  ) +
    plot_layout(heights = c(10, 1))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
