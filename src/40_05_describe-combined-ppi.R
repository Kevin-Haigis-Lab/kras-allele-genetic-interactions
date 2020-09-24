# Descriptive plots of the PPI graph made by combining STRING, HINT, and
# BioPlex2.

GRAPHS_DIR <- "40_05_describe-combined-ppi"
reset_graph_directory(GRAPHS_DIR)

TABLES_DIR <- GRAPHS_DIR
reset_table_directory(TABLES_DIR)


#### ---- Calculate summary statistics ---- ####

component_summary_stats <- combined_ppi_gr %N>%
  convert(to_simple, .select = 1, .clean = TRUE) %E>%
  select(from, to) %N>%
  morph(to_components) %>%
  crystallize() %>%
  mutate(
    num_nodes = map_dbl(graph, igraph::vcount),
    num_edges = map_dbl(graph, igraph::ecount),
    node_degrees = purrr::map(graph, ~ igraph::degree(.x)),
    avg_degree = map_dbl(node_degrees, ~ mean(.x)),
    diameter = map_dbl(graph, ~ igraph::diameter(.x))
  )


#### ---- Plots: components ---- ####

# num nodes vs. avg degree
components_order_size <- component_summary_stats %>%
  ggplot(
    aes(x = log10(num_nodes), y = log10(num_edges + 1))
  ) +
  geom_point(
    aes(size = avg_degree, color = diameter),
    alpha = 0.7
  ) +
  scale_color_viridis_c(option = "D") +
  scale_size_continuous(range = c(0.5, 2)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "graph order (num. nodes)",
    y = "graph size (num. edges)",
    title = "Components of the combined PPIN"
  )
ggsave_wrapper(
  components_order_size,
  plot_path(GRAPHS_DIR, "components_order_size.svg"),
  "small"
)


#### ---- Plots: giant component ---- ####


giant_combined_gr <- combined_ppi_gr %N>%
  jhcutils::get_giant_component()

degree_distribution <- giant_combined_gr %N>%
  mutate(deg = centrality_degree(mode = "all")) %>%
  as_tibble() %>%
  ggplot(aes(x = log2(deg))) +
  geom_density(color = "black", fill = "grey50", alpha = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    x = "log2( node degree )",
    y = "density",
    title = "Degree distribution",
    subtitle = "Distribution of node degree in the giant component of the combined graph"
  )
ggsave_wrapper(
  degree_distribution,
  plot_path(GRAPHS_DIR, "degree_distribution.svg"),
  "small"
)


#### ---- Table: giant component summary stats ---- ####

hubs <- giant_combined_gr %N>%
  convert(to_simple) %>%
  mutate(deg = centrality_degree(mode = "all")) %>%
  as_tibble() %>%
  top_n(n = 10, wt = deg) %>%
  pull(name) %>%
  paste(collapse = ", ")


tibble::tribble(
  ~parameter, ~value,
  "number of nodes", igraph::vcount(giant_combined_gr),
  "numer of edges", igraph::ecount(giant_combined_gr),
  "average degree", mean(igraph::degree(giant_combined_gr)),
  "std. dev. degree", sd(igraph::degree(giant_combined_gr)),
  "diameter", igraph::diameter(giant_combined_gr),
  "hubs", hubs
) %T>%
  write_tsv(
    table_path(TABLES_DIR, "giant_component_stats.tsv")
  )
