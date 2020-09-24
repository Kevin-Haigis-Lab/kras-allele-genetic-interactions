# Compare the subnetworks from LUAD alleles.

source(file.path("src", "40_15_comparing-COAD-allele-subnetworks.R"))

set.seed(0)

GRAPHS_DIR <- "40_17_comparing-PAAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)

PAAD_GRAPH_ANNOTATIONS <- list(
  "ECM" = c("DMD", "COL5A1", "COL18A1", "LAMA1"),
  "calcium signaling" = c(
    "CACNA1C", "CACNA1H", "CALM3", "RYR2",
    "RYR1", "RYR3"
  ),
  "TGFβ signaling" = c("ACVR1B", "ACVR2A", "TGFBR2", "SMAD4", "FKBP1A")
) %>%
  enframe(name = "name", value = "nodes") %>%
  mutate(fill = "grey50")


genes_to_ignore <- c("TTN")
special_nodes <- c(
  "KRAS", "BRAF", "NRAS", "PIK3CA", "TP53",
  "EGFR", "SMAD4", "SMARCA4", "STK11"
)

fxn_grp_file <- table_path(GRAPHS_DIR, "paad_enrichr-output.tsv")

paad_overlap_comparison_plot <- tibble(
  cancer = c("PAAD", "PAAD", "PAAD"),
  allele = c("G12D", "G12R", "G12V")
) %>%
  mutate(
    ppi = map2(cancer, allele,
      get_overlapped_gr,
      min_comp_size = 4,
      ignore_genes = genes_to_ignore
    ),
    data = map(ppi, ~ .x$data),
    ppi = map(ppi, ~ .x$graph),
    # ppi = map(ppi, isolate_kras_subnetwork)
  ) %>%
  make_overlap_comparison_graph() %>%
  get_giant_component() %>%
  print_functional_groups(fxn_grp_file) %>%
  prepare_overlap_comparison_graph_plot(
    special_labels = special_nodes,
    annotation_tib = PAAD_GRAPH_ANNOTATIONS
  ) %>%
  plot_overlap_comparison_graph(
    annotation_tib = PAAD_GRAPH_ANNOTATIONS,
    graph_layout = "stress",
    node_size = 2,
    edge_width = 0.7
  )

ggsave_wrapper(
  paad_overlap_comparison_plot,
  plot_path(GRAPHS_DIR, "paad_overlap_comparison_plot.svg"),
  "large"
)

saveFigRds(paad_overlap_comparison_plot, "paad_overlap_comparison_plot")
