# Compare the subnetworks from LUAD alleles.

source(file.path("src", "40_15_comparing-COAD-allele-subnetworks.R"))

GRAPHS_DIR <- "40_16_comparing-LUAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


isolate_kras_subnetwork <- function(gr) {
    gr %N>%
        morph(to_components) %>%
        mutate(has_kras = any(name == "KRAS")) %>%
        unmorph() %>%
        filter(has_kras)
}



aprior_genes <- wide_genetic_interaction_df %>%
    filter(
        hugo_symbol == "KRAS" | !is.na(KEGG) | !is.na(CGC) | !is.na(BioID)
    ) %>%
    pull(hugo_symbol) %>%
    unique()

LUAD_GRAPH_ANNOTATIONS <- list(
        "WNT signaling" = c("WNT11", "WNT9A", "LRP6", "WNT2", "CTNNB1",
                            "NOTCH1", "FRAT2", "PTPRJ", "TCF7L1"),
        "ARF6" = c("ARF6", "KALRN")
    ) %>%
    enframe(name = "name", value = "nodes") %>%
    mutate(fill = c("grey50", "grey50"))


genes_to_ignore <- c("TTN")
special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53",
                   "EGFR", "SMAD4", "SMARCA4", "STK11")


print_functional_groups <- function(gr, file_name, ignore_genes = NULL) {

    datasources <- c("KEGG_2019_Human", "BioCarta_2016",
                     "GO_Biological_Process_2018", "KEGG_2019_Human",
                     "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human")
    genes <- unique(unlist(igraph::V(gr)$name))
    genes <- genes[!(genes %in% ignore_genes)]
    enrichr_wrapper(genes) %>%
        select(datasource, term, adjusted_p_value, odds_ratio, genes) %>%
        filter(datasource %in% !!datasources) %>%
        filter(adjusted_p_value < 0.05 & odds_ratio > 1.5) %>%
        dplyr::rename(adj_p_value = adjusted_p_value) %>%
        write_tsv(file_name)
    invisible(gr)
}

luad_overlap_comparison_plot <- tibble(
        cancer = c("LUAD", "LUAD"),
        allele = c("G12C", "G12V")
    ) %>%
    mutate(
        ppi = map2(cancer, allele,
                   get_overlapped_gr,
                   min_comp_size = 4,
                   ignore_genes = genes_to_ignore),
        data = map(ppi, ~ .x$data),
        ppi = map(ppi, ~.x$graph),
        ppi = map(ppi, isolate_kras_subnetwork),
        # ppi = map(ppi, ~ .x %N>% filter(name %in% !!aprior_genes)),
        ppi = map(ppi, get_giant_component)
    ) %>%
    make_overlap_comparison_graph() %T>%
    print_functional_groups(
        file_name = table_path(GRAPHS_DIR, "enrichr-output.tsv"),
        ignore_genes = special_nodes
    ) %>%
    plot_overlap_comparison_graph2(
        special_labels = special_nodes,
        annotation_tib = LUAD_GRAPH_ANNOTATIONS,
        node_label_size = 0.5,
        node_label_repel = FALSE,
        node_size = 1
    )

ggsave_wrapper(
    luad_overlap_comparison_plot,
    plot_path(GRAPHS_DIR, "luad_overlap_comparison_plot.svg"),
    "large"
)

saveFigRds(luad_overlap_comparison_plot, "luad_overlap_comparison_plot")
