# Large weakly connected subnetworks from the clusters of the genes.

library(ggraph)

set.seed(0)

ggraph_of_component <- function(gr, idx, cancer = "", cluster_num = "", ...) {

    if (igraph::vcount(gr) == 0) { return() }

    p_title <- glue("{cancer}, cluster {cluster_num} (#{idx})")
    p <- gr %>%
        ggraph(layout = "stress") +
        geom_edge_link(color = "grey25", alpha = 0.5, width = 0.1) +
        geom_node_point(color = "black", size = 0.3) +
        geom_node_text(aes(label = name), family = "arial", repel = TRUE) +
        theme_graph() +
        theme(
            text = element_text("arial")
        ) +
        labs(title = p_title)

    save_path <- plot_path("10_15_linear-modeling-syn-let_ppi-subnetworks",
                           glue("{cancer}_cluster-{cluster_num}_component-{idx}.svg"))
    ggsave_wrapper(p, save_path, "medium")
}
ggraph_of_component <- memoise::memoise(ggraph_of_component)

weakly_connected_components <- function(cancer, gene_cls, hugo_symbols, min_comp_size = 3) {
    hugo_symbols <- unlist(hugo_symbols)
    string_gr %E>%
        filter(combined_score > 900) %N>%
        filter(name %in% !!hugo_symbols) %N>%
        filter(centrality_degree(mode = "all") > 0) %>%
        filter_component_size(min_size = min_comp_size) %>%
        to_components(type = "weak") %>%
        purrr::iwalk(ggraph_of_component, cancer = cancer, cluster_num = gene_cls)
}


depmap_gene_clusters %>%
    group_by(cancer, gene_cls) %>%
    summarise(hugo_symbols = list(hugo_symbol)) %>%
    ungroup() %>%
    purrr::pwalk(weakly_connected_components, min_comp_size = 4)



# depmap_gene_clusters %>%
#     group_by(cancer) %>%
#     summarise(hugo_symbols = list(hugo_symbol)) %>%
#     ungroup() %>%
#     purrr::pwalk(weakly_connected_components, gene_cls = 1, min_comp_size = 5)
