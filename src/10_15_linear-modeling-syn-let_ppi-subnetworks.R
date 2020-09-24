# Large weakly connected subnetworks from the clusters of the genes.

GRAPHS_DIR <- "10_15_linear-modeling-syn-let_ppi-subnetworks"
reset_graph_directory(GRAPHS_DIR)

set.seed(0)


cancers_to_save_for_figures <- c()

save_graph_proto <- function(gg_obj, save_path, cancer) {
  if (cancer %in% cancers_to_save_for_figures) {
    saveFigRds(gg_obj, basename(save_path))
  }
}


ggraph_of_component <- function(gr, idx, cancer = "", cluster_num = "", ...) {
  if (igraph::vcount(gr) == 0) {
    return()
  }

  p_title <- glue("{cancer}, cluster {cluster_num} (#{idx})")
  p <- gr %E>%
    mutate(
      edge_width = centrality_edge_betweenness(directed = FALSE),
      edge_width = scales::rescale(edge_width, to = c(0.6, 2))
    ) %N>%
    mutate(
      node_size = log(centrality_pagerank(directed = FALSE)),
      node_size = scales::rescale(node_size, to = c(1, 3))
    ) %>%
    ggraph(layout = "stress") +
    geom_edge_link(aes(width = edge_width),
      color = "grey40", alpha = 0.5
    ) +
    geom_node_point(aes(size = node_size), color = "grey20") +
    geom_node_text(aes(label = name),
      family = "Arial", size = 2, repel = TRUE
    ) +
    scale_size_identity() +
    scale_edge_width_identity() +
    theme_graph() +
    theme(
      text = element_text("Arial")
    ) +
    labs(title = p_title)

  save_path <- plot_path(
    GRAPHS_DIR,
    glue("{cancer}_cluster-{cluster_num}_component-{idx}.svg")
  )
  ggsave_wrapper(p, save_path, "small")
  save_graph_proto(p, save_path, cancer)
}
ggraph_of_component <- memoise::memoise(ggraph_of_component)


weakly_connected_components <- function(cancer, gene_cls, hugo_symbols,
                                        min_comp_size = 3) {
  hugo_symbols <- unlist(hugo_symbols)
  string_gr %N>%
    filter(name %in% !!hugo_symbols) %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    filter_component_size(min_size = min_comp_size) %>%
    to_components(type = "weak") %>%
    purrr::iwalk(ggraph_of_component,
      cancer = cancer,
      cluster_num = gene_cls
    )
}

depmap_gene_clusters %>%
  group_by(cancer, gene_cls) %>%
  summarise(hugo_symbols = list(hugo_symbol)) %>%
  ungroup() %>%
  purrr::pwalk(weakly_connected_components, min_comp_size = 3)
