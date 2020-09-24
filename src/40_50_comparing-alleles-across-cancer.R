# Comparing the results found for each allele across cancer.

GRAPHS_DIR <- "40_50_comparing-alleles-across-cancer"
reset_graph_directory(GRAPHS_DIR)


#### ---- Prepare `genetic_dependency_df` ---- ####
# A data frame with the allele-to-allele comparison results.

# Get a data frame with two columns `allele1/2_gene_effect` for each allele.
# It is used to merge into the data frame with pairwise results in
# `parse_pariwise_results()`.
gene_effect_average <- function(data, avg_fxn = mean) {
  df1 <- data %>%
    group_by(allele) %>%
    summarise(allele1_gene_effect = avg_fxn(gene_effect)) %>%
    ungroup() %>%
    mutate(allele2_gene_effect = allele1_gene_effect)
}


# Parse the pairwise t-test results of the genetic dependency analysis.
parse_pariwise_results <- function(data, pw_fit, pval_cutoff = Inf) {
  data$allele <- as.character(data$allele)
  gene_effect_avg_df <- gene_effect_average(data)
  ggpubr::compare_means(gene_effect ~ allele,
    data = data,
    method = "t.test", p.adjust.method = "BH"
  ) %>%
    janitor::clean_names() %>%
    filter(p_adj < pval_cutoff) %>%
    dplyr::rename(
      allele1 = group1,
      allele2 = group2,
      pairwise_p_adj = p_adj
    ) %>%
    select(allele1, allele2, pairwise_p_adj) %>%
    left_join(
      select(gene_effect_avg_df, -allele2_gene_effect),
      by = c("allele1" = "allele")
    ) %>%
    left_join(
      select(gene_effect_avg_df, -allele1_gene_effect),
      by = c("allele2" = "allele")
    ) %>%
    mutate(diff_gene_effect = allele1_gene_effect - allele2_gene_effect)
}


# TODO: put this into a ProjectTemplate::cache() when satisfied with it.
genetic_dependency_df <- model1_tib %>%
  filter(rna_pvalue > 0.01) %>%
  mutate(
    aov_p_val = map_dbl(allele_aov, ~ tidy(.x)$p.value[[1]]),
  ) %>%
  filter(aov_p_val < 0.01) %>%
  mutate(
    pariwise_results = map2(data, allele_pairwise, parse_pariwise_results)
  ) %>%
  select(cancer, hugo_symbol, aov_p_val, pariwise_results) %>%
  unnest(pariwise_results)



# Add a duplicate of `df` with the alleles and their values flipped.
add_flipped_alleles <- function(df) {
  df_2 <- df

  df_2$allele1 <- df$allele2
  df_2$allele2 <- df$allele1
  df_2$allele2_gene_effect <- df$allele1_gene_effect
  df_2$allele1_gene_effect <- df$allele2_gene_effect
  df_2$diff_gene_effect <- -1 * df$diff_gene_effect

  bind_rows(df, df_2)
}


comp_genetic_dependency_df <- genetic_dependency_df %>%
  filter(pairwise_p_adj < 0.05) %>%
  add_flipped_alleles() %>%
  filter(allele1 != "WT") %>%
  dplyr::rename(
    allele = allele1,
    comparison_allele = allele2,
    allele_gene_effect = allele1_gene_effect,
    comparison_allele_gene_effect = allele2_gene_effect
  )


#### ---- Merge comutation and genetic dependency ---- ####

comutation_dependency_df <- full_join(
  {
    genetic_interaction_df %>%
      select(cancer, hugo_symbol, allele, p_val, genetic_interaction) %>%
      dplyr::rename(comutation_pval = p_val)
  },
  comp_genetic_dependency_df,
  by = c(
    "cancer" = "cancer",
    "hugo_symbol" = "hugo_symbol",
    "allele" = "allele"
  )
)


comutation_dependency_df %>%
  filter(!is.na(comutation_pval) & !is.na(aov_p_val)) %>%
  filter(hugo_symbol == "WDFY3") %>%
  glimpse()

genetic_interaction_df %>%
  group_by(hugo_symbol, allele) %>%
  filter(n_distinct(cancer) > 1 & n_distinct(genetic_interaction) > 1) %>%
  ungroup() %>%
  arrange(hugo_symbol, allele, cancer) %>%
  select(hugo_symbol, allele, cancer, p_val, genetic_interaction) %>%
  print(n = Inf)

p <- model1_tib %>%
  filter(hugo_symbol == "SETD2") %>%
  select(cancer, data) %>%
  unnest(data) %>%
  ggplot(aes(x = allele, y = gene_effect)) +
  facet_grid(~cancer) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(color = cancer), size = 1) +
  theme_bw(base_size = 7, base_family = "arial")
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "SETD2_all-cancer_box.svg"),
  "medium"
)



assess_damage_predictions <- function(df, num_predictions = 1) {
  df %>%
    mutate(
      num_damaging_preds = is_sift_pred + is_polyphen2_hdiv_pred +
        is_polyphen2_hvar_pred + is_lrt_pred + is_mutation_taster_pred +
        is_mutation_assessor_pred + is_fathmm_pred + is_provean_pred +
        is_fathmm_mkl_coding_pred + is_meta_svm_pred + is_meta_lr_pred +
        is_icgc + is_cosmic + is_clinsig
    ) %>%
    mutate(is_probably_damaging = num_damaging_preds >= num_predictions) %>%
    select(-num_damaging_preds)
}

cancer_coding_av_muts_df %>%
  filter(hugo_symbol == "ABCC9" & cancer == "LUAD") %>%
  assess_damage_predictions() %>%
  filter(is_probably_damaging) %>%
  pull(ras_allele) %>%
  table()
