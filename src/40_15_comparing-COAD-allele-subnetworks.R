# Compare the subnetworks from COAD alleles.

GRAPHS_DIR <- "40_15_comparing-COAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()


#### ---- Comparing of COAD allele overlap networks ---- ####

annotate_merged_graph <- function(gr, df) {
    for (i in 1:nrow(df)) {
        allele_i <- df$allele[[i]]
        nodes <- igraph::V(df$ppi[[i]])$name
        gr <- gr %N>%
            mutate(allele = case_when(
                name %in% !!nodes & allele == "" ~ !!allele_i,
                name %in% !!nodes ~ paste(allele, !!allele_i, sep = ", "),
                TRUE ~ allele
            ))
    }
    return(gr)
}


make_gray_pal <- function(max_ct, start = 0.8, end = 0.5) {
    pal <- grDevices::gray.colors(max_ct - 1, start = start, end = end)
    names(pal) <- paste(seq(2, max_ct), "alleles", sep = " ")
    return(pal)
}


factor_multiple_node_colors <- function(nc) {
    alleles <- unique(nc[!str_detect(nc, "allele")])
    grays <- unique(nc[str_detect(nc, "allele")])
    levels <- c(sort(alleles), sort(grays))
    return(factor(nc, levels = levels))
}


get_allele_count_label <- function(allele) {
    paste(str_count(allele, ",") + 1, "alleles", sep = " ")
}


plot_overlap_comparison_graph <- function(gr, special_labels = NULL) {
    max_ct <- max(str_count(igraph::V(gr)$allele, ",")) + 1
    pal <- c(short_allele_pal,
             make_gray_pal(max_ct),
             "special" = "indianred1")
    p <- gr %N>%
        mutate(
            node_color = case_when(
                name %in% !!special_labels ~ "special",
                str_detect(allele, ",") ~ get_allele_count_label(allele),
                TRUE ~ allele
            ),
            node_color = factor_multiple_node_colors(node_color)
        ) %>%
        ggraph(layout = "kk") +
        geom_edge_link(width = 0.3, color = "gray30", alpha = 0.4) +
        geom_node_point(aes(color = node_color), size = 3) +
        geom_node_text(aes(label = name), size = 2, family = "Arial", repel = TRUE) +
        scale_color_manual(values = pal) +
        theme_graph()
    return(p)
}


make_overlap_comparison_graph <- function(df) {
    merged_gr <- recursive_graph_join(df$ppi) %>%
        convert(to_simple, .clean = TRUE) %>%
        convert(to_undirected, .clean = TRUE) %N>%
        select(name) %E>%
        select(from, to) %N>%
        mutate(allele = "") %>%
        annotate_merged_graph(df)
    return(merged_gr)
}


special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53")

coad_overlap_comparison_plot <- tibble(cancer = c("COAD", "COAD", "COAD"),
       allele = c("G12D", "G12V", "G13D")) %>%
    mutate(
        ppi = purrr::map2(cancer, allele,
                          get_overlapped_gr,
                          min_comp_size = 4, ignore_genes = genes_to_ignore),
        data = purrr::map(ppi, ~ .x$data),
        ppi = purrr::map(ppi, ~.x$graph)
    ) %>%
    make_overlap_comparison_graph() %>%
    plot_overlap_comparison_graph(special_labels = special_nodes)
ggsave_wrapper(
    coad_overlap_comparison_plot,
    plot_path(GRAPHS_DIR, "coad_overlap_comparison_plot.svg"),
    "large"
)
