# Figure 052. #> Comaprison of comutation interactions identified from an
# allele-specific and non-allele-specific analysis.

FIGNUM <- 52

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig52 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_comutation() %+replace%
    theme(
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 6),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      )
    )
}

theme_fig52_venn <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_void(base_size = 7, base_family = "Arial") %+replace%
    theme(
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = tag_margin
      ),
      plot.title = element_text(hjust = 0.5),
      legend.key.size = unit(2, "mm"),
      legend.position = "bottom"
    )
}

theme_fig52_bar <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
  theme_fig52(tag_margin) %+replace%
    theme(
      axis.title.x = element_markdown(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 7,
        family = "Arial",
        hjust = 0.5,
        vjust = 1
      ),
      legend.key.size = unit(2, "mm"),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        vjust = 1
      ),
      legend.box.margin = margin(0, -3, 0, -2, "mm")
    )
}


#### ---- Assemble Venn diagram patch ---- ####

assemble_venn_diagram_patch <- function(cancer,
                                        tag,
                                        tag_margin = margin(-1, -1, -1, -1, "mm")) {
  rc_panel <- read_fig_proto(
    as.character(glue("rc-venn-diagram_{cancer}"))
  ) +
    theme_fig52_venn(tag_margin = tag_margin) +
    labs(tag = tag)
  fish_panel <- read_fig_proto(
    as.character(glue("fish-venn-diagram_{cancer}"))
  ) +
    theme_fig52_venn() %+replace%
    theme(legend.position = "none")

  panel <- ((rc_panel | fish_panel) / guide_area()) +
    plot_layout(
      guides = "collect",
      heights = c(4, 1)
    )
  return(panel)
}


assemble_barplot <- function(cancer, tag, x_label_angle = 0) {
  read_fig_proto(
    as.character(glue("comutation-comparison_bar_{cancer}"))
  ) +
    theme_fig52_bar() +
    labs(tag = tag, x = "*KRAS* alleles")
}


#### ---- A. COAD Venn diagram ---- ####
# A Venn diagram comparing the comutation interactions found when conducting
# an allele specific comutation analysis vs. a non-allele-specific analysis.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_A <- assemble_venn_diagram_patch(
  "COAD",
  tag = "a",
  tag_margin = margin(-1, -1, -1, -1, "mm")
)


#### ---- B. COAD bar plots ---- ####
# A bar plot showing the number of new interactions discovered in the allele-
# specific analysis per KRAS allele.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_B <- assemble_barplot("COAD", "b", x_label_angle = 35)


#### ---- C. LUAD Venn diagram ---- ####
# A Venn diagram comparing the comutation interactions found when conducting
# an allele specific comutation analysis vs. a non-allele-specific analysis.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_C <- assemble_venn_diagram_patch(
  "LUAD",
  tag = "c",
  tag_margin = margin(-1, -1, -1, -1, "mm")
)


#### ---- D. LUAD bar plots ---- ####
# A bar plot showing the number of new interactions discovered in the allele-
# specific analysis per KRAS allele.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_D <- assemble_barplot("LUAD", "d")


#### ---- E. MM Venn diagram ---- ####
# A Venn diagram comparing the comutation interactions found when conducting
# an allele specific comutation analysis vs. a non-allele-specific analysis.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_E <- assemble_venn_diagram_patch(
  "MM",
  tag = "e",
  tag_margin = margin(-1, -1, -1, -1, "mm")
)


#### ---- F. MM bar plots ---- ####
# A bar plot showing the number of new interactions discovered in the allele-
# specific analysis per KRAS allele.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_F <- assemble_barplot("MM", "f")


#### ---- G. PAAD Venn diagram ---- ####
# A Venn diagram comparing the comutation interactions found when conducting
# an allele specific comutation analysis vs. a non-allele-specific analysis.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_G <- assemble_venn_diagram_patch(
  "PAAD",
  tag = "g",
  tag_margin = margin(-1, -1, -1, -1, "mm")
)


#### ---- H. PAAD bar plots ---- ####
# A bar plot showing the number of new interactions discovered in the allele-
# specific analysis per KRAS allele.
# original script: "src/20_67_non-allele-specific-comutation-comparison.R"

panel_H <- assemble_barplot("PAAD", "h")


#### ---- Figure assembly ---- ####

make_cancer_label <- function(cancer) {
  cancer_lbl <- tibble() %>%
    ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = cancer,
      family = "Arial",
      size = 3,
      fontface = "bold"
    ) +
    theme_void(base_family = "Arial", base_size = 7)
  return(cancer_lbl)
}


assemble_cancer_row <- function(venn_panel, bar_panel, cancer, heights = NULL) {
  row_panels <- (
    make_cancer_label(cancer) /
      ((venn_panel | bar_panel) + plot_layout(widths = c(2, 3)))
  ) +
    plot_layout(heights = heights)
  return(row_panels)
}


{
  # COMPLETE FIGURE
  full_figure <- (
    assemble_cancer_row(panel_A, panel_B, "COAD") /
      assemble_cancer_row(panel_C, panel_D, "LUAD", heights = c(1, 10)) /
      assemble_cancer_row(panel_E, panel_F, "MM", heights = c(1, 10)) /
      assemble_cancer_row(panel_G, panel_H, "PAAD", heights = c(1, 10))
  ) +
    plot_layout(
      heights = c(1, 9, 11, 11, 11)
    )

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
