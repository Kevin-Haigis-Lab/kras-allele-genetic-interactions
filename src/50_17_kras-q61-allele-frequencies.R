# Observed vs. predicted KRAS Q61 allele frequency using the trinucleotide
# context.

set.seed(0)

Q61_CODON <- "CAA"

GRAPHS_DIR <- "50_17_kras-q61-allele-frequencies"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


mut_sig_descriptions <- get_mut_sig_descriptions()
artifact_signatures <- get_artifact_signatures()
msg_sigs <- paste0(artifact_signatures, collapse = ", ")
message(glue("The following are artifact signatures: {msg_sigs}"))

# Possible mutations to Q61 codon (ignoring silent and nonsense muts).
possible_q61_alleles <- enumerate_sbs_mutations(Q61_CODON, silent = TRUE) %>%
  filter(amino_acid != "Q") %>%
  filter(amino_acid != "*") %>%
  mutate(kras_allele = as.character(glue("Q61{amino_acid}")))


observed_kras_q61_mutations <- trinucleotide_mutations_df %>%
  filter(hugo_symbol == "KRAS") %>%
  filter(kras_allele != "WT") %>%
  filter(str_detect(kras_allele, "Q61"))

# Assert that all possible mutations were found in the dataset. Can then use the
# mutation's trinucleotide context for further analyses.
stopifnot(
  all(
    possible_q61_alleles$kras_allele %in%
      observed_kras_q61_mutations$amino_acid_change
  )
)

q61_tricontext_mutations <- possible_q61_alleles %>%
  left_join(
    observed_kras_q61_mutations %>%
      select(kras_allele = amino_acid_change, context, tricontext) %>%
      distinct(),
    by = "kras_allele"
  ) %>%
  select(kras_allele, amino_acid, context, tricontext) %>%
  add_column(kras_codon = 61) %>%
  arrange(kras_allele)

knitr::kable(q61_tricontext_mutations)
# |kras_allele |amino_acid |context |tricontext |
# |:-----------|:----------|:-------|:----------|
# |Q61E        |E          |TCA     |T[C>G]A    |
# |Q61H        |H          |CTT     |C[T>A]T    |
# |Q61H        |H          |CTT     |C[T>G]T    |
# |Q61H        |H          |CTT     |C[T>A]T    |
# |Q61H        |H          |CTT     |C[T>G]T    |
# |Q61K        |K          |TCA     |T[C>A]A    |
# |Q61L        |L          |TTG     |T[T>A]G    |
# |Q61P        |P          |TTG     |T[T>G]G    |
# |Q61R        |R          |TTG     |T[T>C]G    |


# Real frequency of Q61 alleles in total populations, just of KRAS mutations,
# and just of KRAS Q61 mutations.

count_cancer_allele_freq <- function(df) {
  df %>%
    distinct(cancer, tumor_sample_barcode, kras_allele) %>%
    count(cancer, kras_allele, name = "n_tsb") %>%
    group_by(cancer) %>%
    mutate(real_allele_freq = n_tsb / sum(n_tsb)) %>%
    ungroup()
}

q61_allele_frequencies_total <- cancer_full_coding_muts_df %>%
  mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
  filter(cancer != "SKCM") %>%
  count_cancer_allele_freq()

q61_allele_frequencies_kras <- cancer_full_coding_muts_df %>%
  mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
  filter(cancer != "SKCM") %>%
  filter(kras_allele != "WT") %>%
  count_cancer_allele_freq()

q61_allele_frequencies_q61 <- cancer_full_coding_muts_df %>%
  mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
  filter(cancer != "SKCM") %>%
  filter(str_detect(kras_allele, "Q61")) %>%
  count_cancer_allele_freq()

# Check frequencies sum to one.
q61_allele_frequencies_q61 %>%
  group_by(cancer) %>%
  summarise(total_freq = sum(real_allele_freq)) %>%
  pull(total_freq) %>%
  purrr::map_lgl(~ abs(.x - 1) < 0.001) %>%
  all() %>%
  stopifnot()


#### ---- Samples with too few mutations ---- ####


MIN_NUM_MUTATION_PER_SAMPLE <- 20

REMOVE_TSB <- trinucleotide_mutations_df %>%
  filter(cancer != "SKCM") %>%
  count(tumor_sample_barcode) %>%
  filter(n < MIN_NUM_MUTATION_PER_SAMPLE) %>%
  pull(tumor_sample_barcode)

message(glue(
  "Removing {length(REMOVE_TSB)} samples because they have too few mutations"
))


#### ---- Samples with all artifact mutational signatures ---- ####


MAXIMUM_ARTIFACT_SIGNATURE <- 0.50
artifact_contribution <- mutational_signatures_df %>%
  filter(cancer != "SKCM") %>%
  filter(!tumor_sample_barcode %in% !!REMOVE_TSB) %>%
  mutate(signature = paste0("sig", signature)) %>%
  total_artifact_contribution(artifact_sigs = artifact_signatures)


ARTIFACT_TSB <- artifact_contribution %>%
  filter(contribution > !!MAXIMUM_ARTIFACT_SIGNATURE) %>%
  pull(tumor_sample_barcode) %>%
  unique()

message(glue(
  "Removing {length(ARTIFACT_TSB)} samples because they have too much Artifact signature."
))

REMOVE_TSB <- c(REMOVE_TSB, ARTIFACT_TSB)

trinucleotide_mutations_df %>%
  filter(cancer != "SKCM") %>%
  filter(!tumor_sample_barcode %in% REMOVE_TSB) %>%
  group_by(cancer, tumor_sample_barcode) %>%
  slice(1) %>%
  ungroup() %>%
  select(cancer, tumor_sample_barcode) %>%
  count(cancer, name = "num_tumor_samples") %>%
  knitr::kable()
# |cancer | num_tumor_samples|
# |:------|-----------------:|
# |COAD   |              1500|
# |LUAD   |               852|
# |MM     |              1185|
# |PAAD   |              1192|


#### ---- Predicted frequencies ---- ####


tricontext_genome_counts <- tricontext_counts_df %>%
  set_names(c("context", "exome", "genome")) %>%
  pivot_longer(-context, names_to = "target", values_to = "genome_count")

q61_tricontext_mutations

q61_cancer_map <- purrr::map_dfr(
  unique(q61_allele_frequencies_q61$cancer),
  ~ q61_tricontext_mutations %>%
    distinct(kras_allele) %>%
    add_column(cancer = .x)
)

cache(
  "kras_q61_predictions",
  depends = c(
    "trinucleotide_mutations_df",
    "q61_tricontext_mutations",
    "q61_allele_frequencies_q61",
    "tricontext_genome_counts",
    "q61_cancer_map",
    "REMOVE_TSB"
  ),
  {
    kras_q61_predictions <- count_tricontext_allele_mutations(
      tricontext_mut_data = trinucleotide_mutations_df,
      allele_df = q61_cancer_map,
      tricontext_counts_df = tricontext_genome_counts,
      kras_tricontexts = q61_tricontext_mutations,
      real_allele_freq = q61_allele_frequencies_q61,
      remove_samples = REMOVE_TSB,
      remove_cancers = c("SKCM")
    ) %>%
      filter_samples_with_no_predictions()
    return(kras_q61_predictions)
  }
)


#### ---- Statistics: boostrapping 95% CI ---- ####


cache("kras_q61_predictions_boot_results",
  depends = c("kras_q61_predictions"),
  {
    set.seed(123)
    kras_q61_predictions_boot_results <- kras_q61_predictions %>%
      group_by(cancer) %>%
      nest() %>%
      mutate(
        boot_res = map2(
          cancer,
          data,
          boot_cancer_expect_frequncies,
          R = 1e3
        ),
        boot_ci = map(boot_res, extract_boot_results)
      )
    return(kras_q61_predictions_boot_results)
  }
)


cancer_q61_expect_frequencies <- kras_q61_predictions %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(data = map(data, calc_expected_frequency)) %>%
  unnest(data)

cancer_q61_expect_frequencies <- kras_q61_predictions_boot_results %>%
  select(cancer, boot_ci) %>%
  unnest(boot_ci) %>%
  right_join(cancer_q61_expect_frequencies, by = c("cancer", "kras_allele"))

# knitr::kable(cancer_q61_expect_frequencies, digits = 3)



#### ---- Check calculations ---- ####

# Check that each TSB sums to 1.
kras_q61_predictions %>%
  group_by(cancer, tumor_sample_barcode) %>%
  check_sum_of_probabilities(prob_col = allele_prob) %>%
  stopifnot()

# Check that each cancer sums to 1.
cancer_q61_expect_frequencies %>%
  group_by(cancer) %>%
  check_sum_of_probabilities(prob_col = expected_allele_frequency) %>%
  stopifnot()


#### ---- Statistics: Chi-Squared ---- ####

q61_chisquared_res <- calc_chisquared_test(
  allele_df = q61_cancer_map,
  expected_freq_df = cancer_q61_expect_frequencies
)


q61_chisquared_res %>%
  filter(adj_p_value > 0.05) %>%
  select(cancer, kras_allele, p_value, adj_p_value) %>%
  knitr::kable(digits = 2)

cancer_q61_expect_frequencies <- cancer_q61_expect_frequencies %>%
  left_join(q61_chisquared_res, by = c("cancer", "kras_allele")) %>%
  select(-parameter, -method) %>%
  rename(
    chi_squared_stat = statistic,
    chi_squared_p_value = p_value,
    chi_squared_adj_p_value = adj_p_value
  )


#### ---- Write tables ---- ####

kras_q61_predictions %>%
  left_join(
    trinucleotide_mutations_df %>%
      distinct(tumor_sample_barcode, cancer, kras_allele) %>%
      rename(observed_kras_allele = kras_allele),
    by = c("tumor_sample_barcode", "cancer")
  ) %>%
  select(-n_tsb) %>%
  mutate_if(is.numeric, scales::label_number(0.001)) %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-q61-predicted-frequencies-per-tumor.tsv"))

cancer_q61_expect_frequencies %>%
  mutate_if(is.numeric, scales::label_number(0.001)) %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-q61-predicted-frequencies.tsv"))


#### ---- Plotting ---- ####

q61_plot_theme <- function() {
  theme_bw(base_size = 7, base_family = "Arial") %+replace%
    theme(
      strip.background = element_blank(),
      axis.ticks = element_blank()
    )
}

plot_cancer_expected_observed_freq <- function(cancer, data) {
  max_val <- max(c(data$expected_allele_frequency, data$observed_allele_frequency))
  data %>%
    ggplot(aes(x = expected_allele_frequency, y = observed_allele_frequency)) +
    geom_hline(yintercept = 0, size = 0.1, color = "grey10") +
    geom_vline(xintercept = 0, size = 0.1, color = "grey10") +
    geom_linerange(aes(xmin = lower_ci, xmax = lower_ci)) +
    geom_point(aes(shape = is_sig), size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    ggrepel::geom_text_repel(
      aes(label = kras_allele),
      size = 2,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    scale_x_continuous(limits = c(0, max_val)) +
    scale_y_continuous(limits = c(0, max_val)) +
    scale_shape_manual(
      values = c(17, 16),
      labels = c("≥ 0.05", "< 0.05"),
      drop = FALSE,
      guide = guide_legend(
        title.position = "top",
        label.position = "top",
        order = 20
      )
    ) +
    coord_equal() +
    q61_plot_theme() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = ifelse(cancer == "PAAD", "right", "none"),
      legend.title = element_markdown(),
      legend.margin = margin()
    ) +
    labs(
      x = "expected",
      y = "observed",
      shape = "χ<sup>2</sup> adj. p-value",
      title = cancer
    )
}

cancer_exp_obs_plots <- cancer_q61_expect_frequencies %>%
  mutate(
    observed_allele_frequency = ifelse(
      is.na(observed_allele_frequency), 0, observed_allele_frequency
    )
  ) %>%
  mutate(is_sig = chi_squared_adj_p_value < 0.05) %>%
  group_by(cancer) %>%
  nest() %>%
  pmap(plot_cancer_expected_observed_freq)

obs_vs_freq_plot <- wrap_plots(cancer_exp_obs_plots, nrow = 1, guides = "collect")
ggsave_wrapper(
  obs_vs_freq_plot,
  plot_path(GRAPHS_DIR, "observed-vs-expected-freq.png"),
  width = 8,
  height = 2.2
)

prediction_dist_plot <- kras_q61_predictions %>%
  ggplot(aes(x = kras_allele, y = allele_prob)) +
  facet_wrap(vars(cancer), nrow = 1) +
  geom_jitter(alpha = 0.2, size = 0.1, width = 0.25) +
  geom_boxplot(outlier.color = NA, alpha = 0.3) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  q61_plot_theme() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown()
  ) +
  labs(x = "*KRAS* allele", y = "probability")

ggsave_wrapper(
  prediction_dist_plot,
  plot_path(GRAPHS_DIR, "prediction-distributions.png"),
  width = 7.5,
  height = 3
)
