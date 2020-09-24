# Fundamental functions for Section 40 in "src" dealing with the integration
# of the comutation and genetic dependency results.


#### ---- Codify the clusters from synthetic lethal analysis ---- ####

# A tibble of the pairwise comparisons between all alleles for the genes in
# the `depmap_gene_clusters`. Set `force = TRUE` to force the remaking of
# of the data frame.
make_depmap_gene_clusters_pairwise_df <- function(force = FALSE) {
  if (exists("depmap_gene_clusters_pairwise_df") & !force) {
    cat("(`depmap_gene_clusters_pairwise_df` already exists)\n")
    return()
  }
  cat("Creating `depmap_gene_clusters_pairwise_df`.\n")

  depmap_gene_clusters_pairwise_df <- depmap_model_workflow_res %>%
    filter_depmap_model_workflow_res() %>%
    inner_join(depmap_gene_clusters, by = c("cancer", "hugo_symbol")) %>%
    select(hugo_symbol, cancer, ova_pairs, gene_cls) %>%
    unnest(ova_pairs) %>%
    select(
      hugo_symbol, cancer, allele, gene_cls,
      p_value, adj_p_value,
      estimate, estimate_allele, estimate_other,
      conf_low, conf_high
    )

  assign("depmap_gene_clusters_pairwise_df",
    depmap_gene_clusters_pairwise_df,
    envir = .GlobalEnv
  )
  ProjectTemplate::cache("depmap_gene_clusters_pairwise_df",
    depends = c("model1_tib", "depmap_gene_clusters")
  )
  invisible()
}


prepare_simple_combined_ppi_gr <- function() {
  if (exists("simple_combined_ppi_gr")) {
    cat("(`simple_combined_ppi_gr` already exists)\n")
    return()
  }

  cat("Creating `simple_combined_ppi_gr`.\n")
  simple_combined_ppi_gr <- convert(combined_ppi_gr, to_simple) %E>%
    mutate(num_source = purrr::map_dbl(
      .orig_data,
      ~ n_distinct(.x$source)
    )) %>%
    select(-.tidygraph_edge_index, -.orig_data) %N>%
    select(-.tidygraph_node_index)

  assign("simple_combined_ppi_gr",
    simple_combined_ppi_gr,
    envir = .GlobalEnv
  )
  ProjectTemplate::cache("simple_combined_ppi_gr",
    depends = c("combined_ppi_gr")
  )
  invisible()
}



# Get the synthetic lethal data for a cancer and allele.
get_synthetic_lethal_data <- function(cancer, allele, adj_p_value = 0.05) {
  depmap_gene_clusters_pairwise_df %>%
    filter(cancer == !!cancer) %>%
    filter(adj_p_value < !!adj_p_value) %>%
    filter(group1 == !!allele | group2 == !!allele) %>%
    select(-cancer)
}


# Get the genetic interaction data for a cancer and allele.
get_genetic_interaction_data <- function(cancer, allele) {
  genetic_interaction_df %>%
    filter(cancer == !!cancer & allele == !!allele) %>%
    select(hugo_symbol, p_val, genetic_interaction) %>%
    mutate(genetic_interaction = ifelse(
      genetic_interaction == "exclusivity",
      "reduced\ncomut.", "increased\ncomut."
    ))
}


# Get a single data frame with all of the dependency and comutation
# data for a cancer and allele.
get_overlapped_df <- function(cancer, allele) {
  dependency_df <- get_synthetic_lethal_data(cancer, allele) %>%
    mutate(
      allele = !!allele,
      other_allele = ifelse(group1 == !!allele, group2, group1),
      comparison = paste(allele, "-", other_allele),
      allele_diff = ifelse(
        group1 == !!allele,
        g1_avg - g1_other_avg,
        g2_avg - g2_other_avg
      )
    ) %>%
    group_by(hugo_symbol) %>%
    slice(1) %>%
    ungroup()
  genetic_df <- get_genetic_interaction_data(cancer, allele)

  df <- full_join(dependency_df, genetic_df, by = "hugo_symbol") %>%
    mutate(interaction_source = case_when(
      is.na(gene_cls) ~ "comutation",
      is.na(genetic_interaction) ~ "dependency",
      TRUE ~ "both"
    ))
  return(df)
}


# Get a PPI annotated with comutation and genetic dependency data.
get_overlapped_gr <- function(cancer, allele, min_comp_size, ignore_genes) {
  merged_df <- get_overlapped_df(cancer, allele)
  gr <- simple_combined_ppi_gr %N>%
    filter(!(name %in% !!ignore_genes)) %>%
    filter(name %in% c(merged_df$hugo_symbol, "KRAS")) %>%
    jhcutils::filter_component_size(min_size = min_comp_size)

  return(list(graph = gr, data = merged_df))
}


print_functional_groups <- function(gr, file_name, ignore_genes = NULL) {
  datasources <- c(
    "KEGG_2019_Human", "BioCarta_2016",
    "GO_Biological_Process_2018", "KEGG_2019_Human",
    "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human"
  )
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



isolate_kras_subnetwork <- function(gr) {
  gr %N>%
    morph(to_components) %>%
    mutate(has_kras = any(name == "KRAS")) %>%
    unmorph() %>%
    filter(has_kras)
}



#### ---- Color palettes ---- ####


kiaa1257_pal <- c(
  "neither" = "grey50",
  "DNAH5" = "#F59237",
  "G12R" = short_allele_pal[["G12R"]],
  "G12R & DNAH5" = "#FF6B72"
)


fkbp1a_edge_types <- c(
  comutation = "dashed",
  dependency = "dotted",
  PPI = "solid"
)

MOD_comut_updown_pal <- comut_updown_pal
names(MOD_comut_updown_pal) <- paste(names(MOD_comut_updown_pal), "comut.")
fkbp1a_edge_pal <- c(
  "PPI regulation" = "grey40",
  "more dep." = synthetic_lethal_pal[["down"]],
  MOD_comut_updown_pal
)

fkbp1a_box_pal <- c(
  short_allele_pal["G12D"],
  short_allele_pal["WT"],
  GPR98 = "#C15953",
  RNF43 = "#479461",
  "GPR98 & RNF43" = "#C7A92E"
)
