# Compare the subnetworks from LUAD alleles.

source(file.path("src", "40_15_comparing-COAD-allele-subnetworks.R"))

GRAPHS_DIR <- "40_17_comparing-PAAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)


genes_to_ignore <- c("TTN")
special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53",
                   "EGFR", "SMAD4", "SMARCA4", "STK11")

paad_overlap_comparison_plot <- tibble(
        cancer = c("PAAD", "PAAD", "PAAD"),
        allele = c("G12D", "G12R", "G12V")
    ) %>%
    mutate(
        ppi = map2(cancer, allele,
                   get_overlapped_gr,
                   min_comp_size = 4,
                   ignore_genes = genes_to_ignore),
        data = map(ppi, ~ .x$data),
        ppi = map(ppi, ~.x$graph),
        ppi = map(ppi, isolate_kras_subnetwork)
    ) %>%
    make_overlap_comparison_graph() %>%
    plot_overlap_comparison_graph2(special_labels = special_nodes)

ggsave_wrapper(
    paad_overlap_comparison_plot,
    plot_path(GRAPHS_DIR, "paad_overlap_comparison_plot.svg"),
    "large"
)


saveFigRds(paad_overlap_comparison_plot, "paad_overlap_comparison_plot")
