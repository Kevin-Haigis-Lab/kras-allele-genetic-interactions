# Figure 054. #> Distributions and statistics for probability of mutational
#  signatures to have caused the various KRAS alleles.

FIGNUM <- 54

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig54 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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

patch_sizes <- list(
  c(2, 5),
  c(1, 10),
  c(1, 4),

  c(1, 10),
  c(1, 9),
  c(1, 10),
  c(1, 10),

  c(1, 10),
  c(1, 10),

  c(1, 6),
  c(1, 10),
  c(1, 10),
  c(1, 10)
)

prepare_ridge_patch <- function(p, idx) {
  p + plot_layout(widths = patch_sizes[[idx]])
}

# source: "src/50_35_mutational-signature-allele-associations.R"
all_panels <- read_fig_proto("ggridge-stats-causation_all-plots") %>%
  imap(prepare_ridge_patch)


#### ---- Figure assembly ---- ####

# styler: off
{
  # COMPLETE FIGURE
  coad_design <- "
  ABC
  "

  luad_design <- "
  AB
  CD
  "

  mm_design <- "
  AB
  "

  paad_design <- "
  AB
  CD
  "

  coad_panels <- wrap_plots(all_panels[1:3], design = coad_design)
  luad_panels <- wrap_plots(all_panels[4:7], design = luad_design)
  mm_panels <- wrap_plots(all_panels[8:9], design = mm_design)
  paad_panels <- wrap_plots(all_panels[10:13], design = paad_design)

  full_figure <- (
    (
      (make_cancer_label("COAD", tag = "a") | coad_panels) +
        plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("LUAD", tag = "b") | luad_panels) +
        plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("MM", tag = "c") | mm_panels) +
        plot_layout(widths = c(1, 30))
    ) /
    (
      (make_cancer_label("PAAD", tag = "d") | paad_panels) +
        plot_layout(widths = c(1, 30))
    )
  ) +
    plot_layout(heights = c(1, 2, 1, 2))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
# styler: on
