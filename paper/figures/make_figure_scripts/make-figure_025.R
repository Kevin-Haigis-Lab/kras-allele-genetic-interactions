# Figure 025. #> PAAD comutation network

FIGNUM <- 25

# > SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 90


theme_fig25 <- function() {
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

# Special theme for graphs from 'ggraph'.
theme_graph_fig25 <- function() {
  theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
      plot.title = element_text(
        size = 7, face = "bold",
        hjust = 0.1, vjust = 5
      ),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(0, 0, 0, 0, "mm")
      ),
      legend.margin = margin(0, 0, 0, 0, "mm"),
      legend.position = "bottom",
      legend.title = element_text(size = 5, hjust = 0.5),
      legend.text = element_text(size = 5, hjust = 0.5),
      plot.margin = margin(0, 0, 0, 0, "mm")
    )
}


#### ---- A. High-level comutation network ---- ####
# A The high-level overview of the comutation network for the *KRAS* alleles
# in PAAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_PAAD") +
  scale_edge_color_manual(
    values = comut_updown_pal,
    guide = guide_legend(
      title = "comutation",
      title.position = "top",
      keywidth = unit(2, "mm"),
      keyheight = unit(1, "mm"),
      nrow = 2,
      label.position = "top",
      order = 1
    )
  ) +
  # scale_color_manual(
  #     values = short_allele_pal,
  #     na.value = NA,
  #     guide = FALSE
  # ) +
  theme_graph_fig25() +
  theme(
    legend.spacing.x = unit(1, "mm"),
    legend.position = c(0.05, 0.05),
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 7),
    legend.margin = margin(-6, 0, 0, 0, "mm")
  ) +
  labs(tag = "a")


#### ---- B. Labeled comutation network of genes in a priori lists ---- ####
# A The labeled comutation network for the *KRAS* alleles only including genes
# in a priori selected gene sets.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"

panel_B <- read_fig_proto(
  "goi_overlap_genetic_interactions_network_PAAD_allLists"
) +
  scale_size_manual(
    values = c(big = 1.6, small = 1.5),
    guide = FALSE
  ) +
  theme_graph_fig25() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(-2, 0, 2, 0, "mm"),
    legend.title = element_markdown()
  ) +
  labs(
    edge_width = "-*log*(p-value)",
    tag = "b"
  )


#### ---- C. Bar plots of log(OR) of select genes ---- ####
# Some genes have two different comutation interactions with multiple KRAS
# alleles. These are bar plots of the log(OR) of a selection of those genes.
# original script: "src/20_41_disagreeing-interactions_logOR-barplot.R"

panel_C <- read_fig_proto(
  "log-odds-ratio_barplot_PAAD.rds"
) +
  facet_wrap(~hugo_symbol, scales = "free", ncol = 1) +
  theme_fig25() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    tag = "c",
    y = "*log*(OR)"
  )



#### ---- Figure assembly ---- ####

{
  set.seed(1)

  # COMPLETE FIGURE
  full_figure <- wrap_elements(full = (
    wrap_elements(full = panel_A) |
      wrap_elements(full = panel_B) |
      panel_C
  ) +
    plot_layout(widths = c(5, 5, 1)))

  save_figure(
    full_figure,
    figure_num = FIGNUM,
    dim = FIG_DIMENSIONS
  )
}
