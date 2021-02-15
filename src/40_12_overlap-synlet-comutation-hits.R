# Individual genes that have both comutation and genetic dependency interactions
# with an allele in a cancer.

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()

#### ---- Overlap of comutation and genetic dependency analysis ---- ####

get_shared_comutation_dependency <- function(cancer, allele) {
  get_overlapped_df(cancer, allele) %>%
    filter(interaction_source == "both") %>%
    mutate(
      cancer = !!cancer,
      genetic_interaction = as.character(genetic_interaction)
    ) %>%
    select(cancer, allele, everything())
}


comut_dep_overlap_tbl <- depmap_gene_clusters_pairwise_df %>%
  distinct(cancer, allele) %>%
  pmap(get_shared_comutation_dependency) %>%
  bind_rows() %>%
  mutate(
    genetic_interaction = str_replace(
      genetic_interaction,
      "\\ncomutation", " comut."
    ),
    comparison = "diff. dep."
  ) %>%
  select(
    cancer, allele, hugo_symbol,
    comparison, adj_p_value,
    genetic_interaction, p_val
  )

knitr::kable(comut_dep_overlap_tbl)

saveFigRds(comut_dep_overlap_tbl, "comut_dep_overlap_tbl.rds")
