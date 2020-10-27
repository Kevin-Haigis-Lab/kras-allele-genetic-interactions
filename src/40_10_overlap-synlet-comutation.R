# Make a data structure for the overlap of the synthetic lethal and genetic
# interaction analyses.

GRAPHS_DIR <- "40_10_overlap-synlet-comutation"
reset_graph_directory(GRAPHS_DIR)

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()


#### ---- Plotting ---- ####

# Dictionary of node shapes to use by genetic interaction.
network_node_shapes <- list(
  "comutation" = 15,
  "dependency" = 17,
  "both" = 18,
  "other" = 19
)

# Dictionary of node colors based on direction of genetic interaction.
network_node_colors <- list(
  "increased\ncomut." = comut_mutex_pal[["comutation"]],
  "reduced\ncomut." = comut_mutex_pal[["exclusivity"]],
  "reduced dep." = synthetic_lethal_pal[["up"]],
  "increased dep." = synthetic_lethal_pal[["down"]],
  "both" = "#ED6165",
  "other" = "#F44C4C"
) %>%
  lighten(factor = 1.4)

# Make the purple even brighter.
network_node_colors["increased dep."] <- lighten(
  network_node_colors[["increased dep."]],
  factor = 1.4
)


# Plot a graph with annotations for comutation and genetic dependencies.
plot_fancy_overlap_ppin <- function(gr, cancer, allele) {
  print(glue("plotting: {cancer} - {allele}"))
  p <- ggraph(gr, layout = "kk") +
    geom_edge_link(
      aes(
        color = num_source,
        width = num_source
      ),
      alpha = 0.7
    ) +
    scale_edge_width_continuous(
      range = c(0.75, 1.5),
      guide = guide_legend(
        title.position = "top",
        label.position = "top",
        order = 3
      )
    ) +
    scale_edge_color_continuous(
      low = "gray80", high = "gray50",
      guide = FALSE
    ) +
    geom_node_point(
      aes(
        shape = interaction_source,
        color = node_color
      ),
      size = 8
    ) +
    geom_node_text(aes(label = name), size = 3, family = "Arial") +
    scale_shape_manual(
      values = unlist(network_node_shapes),
      guide = guide_legend(
        title.position = "top",
        label.position = "top",
        order = 1
      )
    ) +
    scale_color_manual(
      values = network_node_colors,
      guide = guide_legend(
        title.position = "top",
        label.position = "top",
        order = 2
      )
    ) +
    theme_graph(base_family = "Arial", base_size = 8) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      title = glue("Interactions on the PPI for {allele} in {cancer}"),
      shape = "type of interaction",
      color = "details of interaction"
    )
  return(p)
}


# Assign the node shape and color based on a combination of the comutation
# and the dependency analysis.
assign_node_color_and_shape <- function(gr) {
  new_gr <- gr %N>%
    mutate(
      genetic_interaction = ifelse(
        is.na(genetic_interaction), NA_character_, genetic_interaction
      ),
      interaction_source = ifelse(
        is.na(interaction_source), "other", interaction_source
      ),
      synlet_direction = ifelse(
        allele_diff < 0, "increased dep.", "reduced dep."
      ),
      node_shape = network_node_shapes[interaction_source],
      node_color = case_when(
        interaction_source == "comutation" ~ genetic_interaction,
        interaction_source == "dependency" ~ synlet_direction,
        interaction_source == "other" ~ "other",
        interaction_source == "both" ~ "both"
      )
    )
  return(new_gr)
}


fancy_overlap_ppin <- function(cancer, allele,
                               min_comp_size = 4,
                               ignore_genes = c()) {
  set.seed(0)

  print(glue("beginning: {cancer} - {allele}"))

  res <- get_overlapped_gr(cancer, allele, min_comp_size, ignore_genes)
  gr <- res[["graph"]]
  df <- res[["data"]]

  if (igraph::vcount(gr) == 0 | igraph::ecount(gr) == 0) {
    return(NULL)
  }

  save_graph_plot <- function(name, gr_plot, ...) {
    print(glue("saving: {cancer} - {allele}"))

    save_path <- plot_path(
      GRAPHS_DIR,
      paste0("overlap_ppi_", cancer, "_", allele, "_", name, ".svg")
    )
    ggsave_wrapper(gr_plot, save_path, width = 10, height = 8)
  }

  gr_components <- gr %N>%
    left_join(df, by = c("name" = "hugo_symbol")) %>%
    assign_node_color_and_shape() %>%
    morph(to_components) %>%
    crystallize() %>%
    mutate(
      gr_plot = purrr::map(graph, plot_fancy_overlap_ppin,
        cancer = !!cancer, allele = !!allele
      )
    ) %>%
    pwalk(save_graph_plot)
}


clustered_fancy_overlap_ppin <- function(cancer, allele,
                                         min_comp_size = 4,
                                         ignore_genes = c()) {
  set.seed(0)
  print(glue("beginning (clustered): {cancer} - {allele}"))

  res <- get_overlapped_gr(cancer, allele, min_comp_size, ignore_genes)
  gr <- res[["graph"]]
  df <- res[["data"]]

  if (igraph::vcount(gr) == 0 | igraph::ecount(gr) == 0) {
    return(NULL)
  }

  save_graph_plot <- function(name, gr_plot, ...) {
    print(glue("saving (clustered): {cancer} - {allele}"))

    save_path <- plot_path(
      GRAPHS_DIR,
      paste0("clustered_ppi_", cancer, "_", allele, "_", name, ".svg")
    )
    ggsave_wrapper(gr_plot, save_path, width = 10, height = 8)
  }

  gr_components <- gr %N>%
    left_join(df, by = c("name" = "hugo_symbol")) %>%
    assign_node_color_and_shape() %N>%
    morph(to_components) %>%
    mutate(cls = group_spinglass()) %>%
    unmorph() %E>%
    filter(.N()$cls[from] == .N()$cls[to]) %N>%
    morph(to_components) %>%
    crystallize() %>%
    mutate(
      gr_plot = purrr::map(graph, plot_fancy_overlap_ppin,
        cancer = !!cancer, allele = !!allele
      )
    ) %>%
    pwalk(save_graph_plot)
}


genes_to_ignore <- c("TTN")

depmap_gene_clusters_pairwise_df %>%
  select(cancer, group1, group2) %>%
  unique() %>%
  group_by(cancer) %>%
  summarise(allele = list(unique(c(unlist(group1), unlist(group2))))) %>%
  ungroup() %>%
  unnest(allele) %>%
  pwalk(fancy_overlap_ppin, ignore_genes = genes_to_ignore) %>%
  pwalk(clustered_fancy_overlap_ppin, ignore_genes = genes_to_ignore)
