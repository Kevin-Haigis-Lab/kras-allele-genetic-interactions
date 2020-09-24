
# Code for generating the information and figures for KH's R01 resubmission.

library(ggraph)

# For reproducibility.
set.seed(0)

# Set-up directory to which images are saved.
SAVE_DIR <- "90_15_kh-resubmission"
reset_graph_directory(SAVE_DIR)

# Parameterize the cancer in case there are more in the future.
CANCER <- "COAD"

# Genes of interest.
goi_genes <- unique(genes_of_interest_df$hugo_symbol)

make_thick_network_plot <- function(gr) {
  p <- gr %>%
    ggraph(layout = "nicely") +
    geom_edge_link(
      aes(color = genetic_interaction),
      width = 0.3
    ) +
    scale_edge_color_manual(values = comut_mutex_pal) +
    geom_node_point(
      aes(
        color = node_color,
        size = node_size
      ),
    ) +
    geom_node_text(
      aes(
        label = node_label,
        size = label_size
      ),
      repel = TRUE,
      family = "Arial"
    ) +
    scale_color_manual(values = short_allele_pal, na.value = NA) +
    scale_size_identity() +
    theme_graph(base_family = "Arial", base_size = 8) +
    theme(
      plot.title = element_blank(),
      legend.position = "none"
    )
  return(p)
}



gr_plot <- genetic_interaction_gr %E>%
  filter(cancer == !!CANCER) %N>%
  filter(centrality_degree(mode = "all") > 0) %>%
  mutate(
    node_label = ifelse(is_kras, name, NA),
    node_label = str_remove_all(node_label, "KRAS_"),
    node_color = node_label,
    node_size = 2,
    label_size = 5
  ) %>%
  make_thick_network_plot() +
  labs(
    title = glue("Genetic interactions in {CANCER}"),
    color = "KRAS allele",
    edge_color = "genetic\ninteraction"
  )
save_path <- plot_path(
  SAVE_DIR,
  glue("genetic_interaction_network_{CANCER}_thick.svg")
)
ggsave_wrapper(gr_plot, save_path, size = "medium")


gr_plot <- genetic_interaction_gr %E>%
  filter(cancer == !!CANCER & genetic_interaction == "comutation") %N>%
  filter(centrality_degree(mode = "all") > 0) %>%
  mutate(
    node_label = ifelse(is_kras | name %in% !!goi_genes, name, NA),
    node_label = str_remove_all(node_label, "KRAS_"),
    node_color = case_when(
      is_kras ~ node_label,
      name %in% !!goi_genes ~ "Other"
    ),
    node_size = case_when(
      is_kras ~ 2,
      name %in% !!goi_genes ~ 0.5
    ),
    label_size = case_when(
      is_kras ~ 5,
      name %in% !!goi_genes ~ 2.5
    )
  ) %>%
  make_thick_network_plot() +
  labs(
    title = glue("Genetic interactions in {CANCER}"),
    color = "KRAS allele",
    edge_color = "genetic\ninteraction"
  )

save_path <- plot_path(
  SAVE_DIR,
  glue("genetic_interaction_network_{CANCER}_thick_comutation_goi.svg")
)
ggsave_wrapper(gr_plot, save_path, size = "medium")


#### ---- Prioritizing genes with damaging mutations ---- ####

# Just the COAD samples.
cancer_coding_av_muts_df_COAD <- filter(cancer_coding_av_muts_df, cancer == !!CANCER)

# For a data frame that is subset of `cancer_coding_av_muts_df`, return just
# the rows that are predicted to be harmful by a predictive algorithm or
# experimental system.
filter_damaging_mutations <- function(df) {
  filter(
    df,
    is_sift_pred | is_polyphen2_hdiv_pred | is_polyphen2_hvar_pred |
      is_lrt_pred | is_mutation_taster_pred | is_mutation_assessor_pred |
      is_fathmm_pred | is_provean_pred | is_fathmm_mkl_coding_pred |
      is_meta_svm_pred | is_meta_lr_pred | is_gerp |
      is_phast_cons20way_mammalian |
      is_si_phy_29way_log_odds | is_icgc | is_cosmic | is_clinsig
  )
}

# For a data frame that is subset of `cancer_coding_av_muts_df`, return the
# number of predictive metrics that predict damage.
count_damaging_predictions <- function(df) {
  df[is.na(df)] <- FALSE
  df <- mutate(df,
    num_predicted =
      is_sift_pred + is_polyphen2_hdiv_pred + is_polyphen2_hvar_pred +
        is_lrt_pred + is_mutation_taster_pred + is_mutation_assessor_pred +
        is_fathmm_pred + is_provean_pred + is_fathmm_mkl_coding_pred +
        is_meta_svm_pred + is_meta_lr_pred + is_gerp +
        is_phast_cons20way_mammalian +
        is_si_phy_29way_log_odds + is_icgc + is_cosmic + is_clinsig
  )
  return(df$num_predicted)
}

# For a gene and (optional) KRAS allele, find the number of mutations predicted
# to be damaging by at least `min_predictions_per_mutation` algorithms.
num_predictions_of_damage <- function(gene,
                                      kras_allele = NULL,
                                      min_predictions_per_mutation = 1) {
  df <- cancer_coding_av_muts_df_COAD %>%
    filter(hugo_symbol == gene)

  if (nrow(df) == 0) {
    return(0)
  }

  if (!is.null(kras_allele)) {
    df <- filter(df, ras_allele %in% !!kras_allele)
  }

  if (nrow(df) == 0) {
    return(0)
  }

  df$num_predicted_damaging <- count_damaging_predictions(df)
  df %<>% filter(num_predicted_damaging >= !!min_predictions_per_mutation)

  return(nrow(df))
}


# Just the total number of mutations found between a gene and
# a (optional) KRAS allele.
number_of_mutations <- function(gene, kras_allele = NULL) {
  df <- cancer_coding_av_muts_df_COAD %>%
    filter(hugo_symbol == gene)

  if (nrow(df) == 0) {
    return(0)
  }

  if (!is.null(kras_allele)) {
    df <- filter(df, ras_allele %in% !!kras_allele)
  }

  return(nrow(df))
}


KRAS_ALLELES_TO_USE <- paste0("KRAS_", c("G12D", "G12V", "G12C", "G13D", "A146T"))
cgc_genes <- unique(genes_of_interest_df$hugo_symbol[genes_of_interest_df$source == "CGC"])

gr_damaging <- genetic_interaction_gr %E>%
  filter(kras_allele %in% !!KRAS_ALLELES_TO_USE) %>%
  filter(cancer == !!CANCER & genetic_interaction == "comutation") %N>%
  filter(!is_kras | name %in% !!KRAS_ALLELES_TO_USE) %>%
  filter(centrality_degree(mode = "all") > 0) %E>%
  mutate(
    num_damaged = purrr::map2_dbl(hugo_symbol, kras_allele,
      num_predictions_of_damage,
      min_predictions_per_mutation = 5
    ),
    num_mutated = purrr::map2_dbl(hugo_symbol, kras_allele, number_of_mutations)
  ) %E>%
  mutate(fraction_damaging = num_damaged / num_mutated) %>%
  filter(fraction_damaging >= 0.70 |
    hugo_symbol %in% cgc_genes) %N>%
  filter(centrality_degree(mode = "all") > 0)

as_tibble(gr_damaging, active = "nodes") %>%
  filter(!is_kras) %>%
  pull(name) %>%
  n_distinct()

# > 42

as_tibble(gr_damaging, active = "edges") %>%
  count(kras_allele)

# > Allele    no. genes
# > A146T             4
# > G12C              7
# > G12D             17
# > G12V              8
# > G13D              8

gr_damaging_plot <- gr_damaging %N>%
  mutate(
    node_label = str_remove_all(name, "KRAS_"),
    node_color = case_when(
      is_kras ~ node_label,
      TRUE ~ "Other"
    ),
    node_size = case_when(
      is_kras ~ 2,
      TRUE ~ 0.1
    ),
    label_size = case_when(
      is_kras ~ 5,
      TRUE ~ 2.5
    )
  ) %>%
  make_thick_network_plot() +
  labs(
    title = glue("Genetic interactions in {CANCER}"),
    color = "KRAS allele",
    edge_color = "genetic\ninteraction"
  )

save_path <- plot_path(
  SAVE_DIR,
  glue("genetic_interaction_network_{CANCER}_thick_comutation_goi_damaging.svg")
)
ggsave_wrapper(gr_damaging_plot, save_path, size = "medium")
