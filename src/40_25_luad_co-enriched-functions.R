# Analysis of two similar functions enriched in the comutation and genetic
# dependency analysis of LUAD G12C samples.

GRAPHS_DIR <- "40_25_luad_co-enriched-functions.R"
reset_graph_directory(GRAPHS_DIR)



#### ---- Enriched functions and genes ---- ####

dependency_enrchfxn <- list(
  geneset_family = "BIOCARTA",
  geneset = "ERK_PATHWAY"
)

comutation_enrchfxn <- list(
  geneset_family = "GO_Biological_Process_2018",
  geneset = "positive regulation of MAPK cascade"
)


comutation_enrchfxn_genes <- enrichr_tib %>%
  filter(cancer == "LUAD" & allele == "G12C") %>%
  select(-gene_list) %>%
  unnest(enrichr_res) %>%
  mutate(term = standardize_enricher_terms(term)) %>%
  filter(
    datasource == comutation_enrchfxn$geneset_family &
      term == comutation_enrchfxn$geneset
  ) %>%
  mutate(overlap_genes = str_split(genes, ";")) %>%
  select(-genes) %>%
  group_by(term) %>%
  mutate(term_genes = list(unique(unlist(overlap_genes)))) %>%
  pull(term_genes) %>%
  unlist() %>%
  unique()


dependency_enrchfxn_genes <- gsea_df %>%
  filter(cancer == "LUAD" & allele == "G12C") %>%
  filter(
    gene_set_family == dependency_enrchfxn$geneset_family &
      gene_set == dependency_enrchfxn$geneset
  ) %>%
  select(cancer, allele, name) %>%
  pmap(get_geneset_enrichment_results) %>%
  bind_rows() %>%
  filter(core_enrichment) %>%
  pull(probe) %>%
  unlist() %>%
  unique()


enrchfxn_genes <- unique(unlist(c(
  comutation_enrchfxn_genes,
  dependency_enrchfxn_genes
)))


#### ---- Reduced KEGG pathways ---- ####


KEGG_PATHWAYS_TO_USE <- c(
  "mapk_signaling_pathway",
  "ras_signaling_pathway",
  "pi3k_akt_signaling_pathway",
  "wnt_signaling_pathway"
)


kegg_graph_join <- function(x, y) {
  x <- x %N>% select(name)
  y <- y %N>% select(name)

  graph_join(x, y, by = "name")
}


for (pwname in KEGG_PATHWAYS_TO_USE) {
  print(pwname)
  kegg_pathway_grs[[pwname]] %N>%
    as_tibble() %>%
    group_by(name) %>%
    summarise(num = n()) %>%
    ungroup() %>%
    filter(num > 1) %>%
    print()
}



kegg_gr <- kegg_pathway_grs[[KEGG_PATHWAYS_TO_USE[[1]]]]
for (pwname in KEGG_PATHWAYS_TO_USE[-1]) {
  kegg_gr <- kegg_graph_join(kegg_gr, kegg_pathway_grs[[pwname]])
}


# Get all pairs of the nodes in `gs`.
get_all_node_pairs <- function(gs) {
  combn(gs, 2) %>%
    t() %>%
    as.data.frame(strings_as_factors = FALSE) %>%
    as_tibble() %>%
    set_colnames(c("from", "to"))
}

# Returns a boolean as to whether a path exists between `from` and `to`.
# `from` and `to` must be node indices.
path_does_exist <- function(gr, from, to) {
  !any(is.infinite(igraph::distances(gr, from, to)))
}


# Get the nodes along the shortest path(s) between `from` and `to`.
# Path the `from` and `to` as node names.
get_nodes_in_shortest_path <- function(from, to, gr) {
  from_idx <- get_node_index(gr, name == from)
  to_idx <- get_node_index(gr, name == to)

  # Skip if there are no shortest paths.
  if (!path_does_exist(gr, from_idx, to_idx)) {
    return(c(from, to))
  }

  # Get names of the nodes on the shortest path.
  names(unlist(igraph::shortest_paths(gr, from, to, mode = "all")$vpath))
}


# Shrink the graph to only include the nodes that lie on the geodesic path
# of the main nodes.
minimize_graph_to_shortest_paths <- function(gr, main_nodes) {
  main_nodes <- main_nodes[main_nodes %in% igraph::V(gr)$name]
  node_pairs <- get_all_node_pairs(main_nodes)
  genes_on_paths <- pmap(node_pairs, get_nodes_in_shortest_path, gr = gr) %>%
    unlist() %>%
    unique()
  mod_gr <- gr %N>%
    filter(name %in% genes_on_paths)
  return(mod_gr)
}
minimize_graph_to_shortest_paths <- memoise::memoise(minimize_graph_to_shortest_paths)


get_reduced_subnetwork <- function(base_gr, core_nodes) {
  gr <- base_gr %>%
    get_giant_component() %>%
    mutate(
      node_ctrlty = centrality_betweenness(directed = FALSE),
      node_ctrlty = ifelse(
        name %in% core_nodes, max(node_ctrlty), node_ctrlty
      )
    )

  if (!all(core_nodes %in% igraph::V(gr)$name)) {
    message("Not all nodes in giant component -> returned original graph")
    return(base_gr)
  }

  MAX_CTRLTY <- max(igraph::V(gr)$node_ctrlty)
  while (!all(igraph::V(gr)$name %in% core_nodes)) {
    mod_gr <- gr %N>% filter(node_ctrlty > min(node_ctrlty))

    if (igraph::count_components(mod_gr) == 1) {
      gr <- mod_gr
    } else {
      gr <- gr %N>% mutate(node_ctrlty = ifelse(
        node_ctrlty == min(node_ctrlty), !!MAX_CTRLTY, node_ctrlty
      ))
    }

    # Stop if all of the nodes have the maximum centrality.
    if (all(igraph::V(gr)$node_ctrlty == MAX_CTRLTY)) break
  }

  return(gr)
}

test_genes <- c("KRAS", "NRAS", "EGFR")
kegg_gr %N>%
  filter(name %in% test_genes)


gr_plot <- kegg_gr %>%
  minimize_graph_to_shortest_paths(enrchfxn_genes) %>%
  get_reduced_subnetwork(enrchfxn_genes) %N>%
  mutate(is_enriched = name %in% enrchfxn_genes) %>%
  ggraph("nicely") +
  geom_edge_link() +
  geom_node_point() +
  geom_node_text(aes(label = name, color = is_enriched), repel = TRUE) +
  theme_graph()
ggsave_wrapper(
  gr_plot,
  plot_path(GRAPHS_DIR, "test_plot.svg"),
  "medium"
)
