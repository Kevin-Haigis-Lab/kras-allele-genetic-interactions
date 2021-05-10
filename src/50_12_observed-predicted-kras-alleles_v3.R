# Observed vs. predicted KRAS allele frequency using the trinucleotide context.
# This is a new version for "src/50_10_observed-predicted-kras-alleles.R"
# because I am not confident that it is bug-free.
#
# Also, GM and I discussed the statistical testing and decided to calculate and
# R-squared value and use the Chi squared test instead of a binomial.

set.seed(0)

library(ggrepel)

GRAPHS_DIR <- "50_12_observed-predicted-kras-alleles_v3"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)

mut_sig_descriptions <- get_mut_sig_descriptions()
artifact_signatures <- get_artifact_signatures()
msg_sigs <- paste0(artifact_signatures, collapse = ", ")
message(glue("The following are artifact signatures: {msg_sigs}"))


oncogenic_alleles <- unique(kras_trinucleotide_contexts$kras_allele)
oncogenic_alleles <- oncogenic_alleles[!oncogenic_alleles %in% c("K117N")]
message("The following are the oncogenic KRAS alleles:")
print(oncogenic_alleles)


#### ---- KRAS allele frequency ---- ####


real_kras_allele_freq <- cancer_full_coding_muts_df %>%
  filter(cancer != "SKCM") %>%
  distinct(cancer, tumor_sample_barcode, ras_allele) %>%
  mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
  group_by(cancer) %>%
  summarise(
    kras_freq = list(
      calc_frequency_of_oncogenic_kras(
        kras_allele,
        oncogenic_alleles = oncogenic_alleles
      )
    )
  ) %>%
  unnest(kras_freq) %>%
  rename(real_allele_freq = frequency)

wegs_kras_allele_ct <- trinucleotide_mutations_df %>%
  distinct(cancer, tumor_sample_barcode, kras_allele) %>%
  count(cancer, kras_allele, name = "num_allele_tsb")



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

artifact_contribution_plt <- artifact_contribution %>%
  filter(contribution > 0) %>%
  ggplot(aes(x = contribution)) +
  geom_histogram(breaks = seq(0, 1, 0.05)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 7, base_family = "Arial") +
  labs(
    x = "contribution of Artifact",
    y = "count"
  )
ggsave_wrapper(
  artifact_contribution_plt,
  plot_path(GRAPHS_DIR, "artifact_contribution_plt.svg"),
  "small"
)

ARTIFACT_TSB <- artifact_contribution %>%
  filter(contribution > !!MAXIMUM_ARTIFACT_SIGNATURE) %>%
  pull(tumor_sample_barcode) %>%
  unique()

message(glue(
  "Removing {length(ARTIFACT_TSB)} samples because they have too much Artifact signature."
))

REMOVE_TSB <- c(REMOVE_TSB, ARTIFACT_TSB)

# Remaining number of samples per cancer.
trinucleotide_mutations_df %>%
  filter(cancer != "SKCM") %>%
  filter(!tumor_sample_barcode %in% REMOVE_TSB) %>%
  group_by(cancer, tumor_sample_barcode) %>%
  slice(1) %>%
  ungroup() %>%
  select(cancer, tumor_sample_barcode) %>%
  count(cancer, name = "num_tumor_samples") %>%
  knitr::kable()

# > |cancer | num_tumor_samples|
# > |:------|-----------------:|
# > |COAD   |              1500|
# > |LUAD   |               852|
# > |MM     |              1185|
# > |PAAD   |              1192|


#### ---- Calculate predicted frequency ---- ####

tricontext_genome_counts <- tricontext_counts_df %>%
  set_names(c("context", "exome", "genome")) %>%
  pivot_longer(-context, names_to = "target", values_to = "genome_count")


# Limiting alleles to only those found in the cancer at a certain frequency.
MINIMUM_ALLELE_FREQ <- 0.03
TOP_ALLELES_PER_CANCER <- 5

cache(
  "alleles_for_each_cancer_obs_vs_pred",
  depends = c(
    "MINIMUM_ALLELE_FREQ",
    "TOP_ALLELES_PER_CANCER",
    "real_kras_allele_freq"
  ),
  {
    alleles_for_each_cancer_obs_vs_pred <- real_kras_allele_freq %>%
      filter(cancer != "SKCM") %>%
      left_join(wegs_kras_allele_ct, by = c("cancer", "kras_allele")) %>%
      arrange(cancer, desc(real_allele_freq)) %>%
      group_by(cancer) %>%
      mutate(cancer_rank = 1:n()) %>%
      ungroup() %>%
      filter(
        real_allele_freq >= MINIMUM_ALLELE_FREQ |
          cancer_rank <= TOP_ALLELES_PER_CANCER
      ) %>%
      filter(num_allele_tsb >= 10) %>%
      select(cancer, kras_allele)
    return(alleles_for_each_cancer_obs_vs_pred)
  }
)

alleles_for_each_cancer_obs_vs_pred %>%
  group_by(cancer) %>%
  summarise(alleles = paste0(kras_allele, collapse = ", ")) %>%
  knitr::kable(format = "simple")
# > cancer   alleles
# > -------  -----------------------------------------------
# > COAD     G12D, G12V, G13D, A146T, G12C, G12A, G12S
# > LUAD     G12C, G12V, G12D, G12A, G13C
# > MM       Q61H, G13D, G12D, G12V, G12A, Q61R, G12R, Q61L
# > PAAD     G12D, G12V, G12R, Q61H, G12C

cache(
  "kras_allele_predictions",
  depends = c(
    "trinucleotide_mutations_df",
    "alleles_for_each_cancer_obs_vs_pred",
    "REMOVE_TSB"
  ),
  {
    kras_allele_predictions <- count_tricontext_allele_mutations(
      tricontext_mut_data = trinucleotide_mutations_df,
      allele_df = alleles_for_each_cancer_obs_vs_pred,
      tricontext_counts_df = tricontext_genome_counts,
      kras_tricontexts = kras_trinucleotide_contexts,
      real_allele_freq = real_kras_allele_freq,
      remove_samples = REMOVE_TSB,
      remove_cancers = c("SKCM")
    ) %>%
      filter_samples_with_no_predictions()
    return(kras_allele_predictions)
  }
)



# Using all oncogenic alleles found in any of the cancers.
cache(
  "all_observed_alleles_obs_vs_pred",
  depends = c(
    "kras_trinucleotide_contexts",
    "alleles_for_each_cancer_obs_vs_pred"
  ),
  {
    all_observed_alleles_obs_vs_pred <- expand.grid(
      kras_allele = sort(unique(kras_trinucleotide_contexts$kras_allele)),
      cancer = sort(unique(alleles_for_each_cancer_obs_vs_pred$cancer)),
      stringsAsFactors = FALSE
    ) %>%
      filter(!str_detect(kras_allele, "117")) %>%
      as_tibble()
    return(all_observed_alleles_obs_vs_pred)
  }
)


cache(
  "all_kras_allele_predictions",
  depends = c(
    "trinucleotide_mutations_df",
    "all_observed_alleles_obs_vs_pred",
    "REMOVE_TSB"
  ),
  {
    all_kras_allele_predictions <- count_tricontext_allele_mutations(
      tricontext_mut_data = trinucleotide_mutations_df,
      allele_df = all_observed_alleles_obs_vs_pred,
      tricontext_counts_df = tricontext_genome_counts,
      kras_tricontexts = kras_trinucleotide_contexts,
      real_allele_freq = real_kras_allele_freq,
      remove_samples = REMOVE_TSB,
      remove_cancers = c("SKCM")
    ) %>%
      filter_samples_with_no_predictions()
    return(all_kras_allele_predictions)
  }
)



#### ---- Statistics: boostrapping 95% CI ---- ####

## For KRAS alleles found in each cancer as a high frequency.

# Results from bootstrapping the samples used for the calculation of the
# predicted KRAS alleles.
cache("kras_allele_predictions_boot_results",
  depends = c("kras_allele_predictions"),
  {
    set.seed(123)
    kras_allele_predictions_boot_results <- kras_allele_predictions %>%
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
    return(kras_allele_predictions_boot_results)
  }
)

# Real values
cancer_expect_frequencies <- kras_allele_predictions %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(data = map(data, calc_expected_frequency)) %>%
  unnest(data)


cancer_expect_frequencies <- kras_allele_predictions_boot_results %>%
  select(cancer, boot_ci) %>%
  unnest(boot_ci) %>%
  right_join(cancer_expect_frequencies, by = c("cancer", "kras_allele"))

knitr::kable(cancer_expect_frequencies, digits = 3)
cancer_expect_frequencies %>%
  mutate_if(is.numeric, scales::label_number(0.001)) %>%
  save_supp_data(6, "pred vs obs KRAS alleles")



## For all KRAS alleles for in any cancer.

# Results from bootstrapping the samples used for the calculation of the
# predicted KRAS alleles.
cache("all_kras_allele_predictions_boot_results",
  depends = c("all_kras_allele_predictions"),
  {
    set.seed(123)
    all_kras_allele_predictions_boot_results <- all_kras_allele_predictions %>%
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
    return(all_kras_allele_predictions_boot_results)
  }
)


all_expect_frequencies <- all_kras_allele_predictions %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(data = map(data, calc_expected_frequency)) %>%
  unnest(data)

all_expect_frequencies <- all_kras_allele_predictions_boot_results %>%
  select(cancer, boot_ci) %>%
  unnest(boot_ci) %>%
  right_join(all_expect_frequencies, by = c("cancer", "kras_allele"))
knitr::kable(all_expect_frequencies, digits = 3)
all_expect_frequencies %>%
  mutate_if(is.numeric, scales::label_number(0.001)) %>%
  save_supp_data(10, "pred vs obs all KRAS alleles")


#### ---- Check calculations ---- ####

# Check that each TSB sums to 1.
kras_allele_predictions %>%
  group_by(cancer, tumor_sample_barcode) %>%
  check_sum_of_probabilities(prob_col = allele_prob) %>%
  stopifnot()
all_kras_allele_predictions %>%
  group_by(cancer, tumor_sample_barcode) %>%
  check_sum_of_probabilities(prob_col = allele_prob) %>%
  stopifnot()

# Check that each cancer sums to 1.
cancer_expect_frequencies %>%
  group_by(cancer) %>%
  check_sum_of_probabilities(prob_col = expected_allele_frequency) %>%
  stopifnot()
all_expect_frequencies %>%
  group_by(cancer) %>%
  check_sum_of_probabilities(prob_col = expected_allele_frequency) %>%
  stopifnot()



#### ---- Statistics: R-squared ---- ####


cancer_rsquareds <- cancer_expect_frequencies %>%
  group_by(cancer) %>%
  summarise(model_fit = list(
    calc_obs_pred_rsquared(
      observed_allele_frequency,
      expected_allele_frequency
    )
  )) %>%
  unnest(model_fit) %>%
  ungroup() %>%
  select(cancer, r_squared, adj_r_squared, p_value) %>%
  rename(model_p_value = p_value)
cancer_correlations <- cancer_expect_frequencies %>%
  group_by(cancer) %>%
  summarise(model_fit = list(
    calc_obs_pred_correlation(
      observed_allele_frequency,
      expected_allele_frequency
    )
  )) %>%
  unnest(model_fit)


all_rsquareds <- all_expect_frequencies %>%
  group_by(cancer) %>%
  summarise(model_fit = list(
    calc_obs_pred_rsquared(
      observed_allele_frequency,
      expected_allele_frequency
    )
  )) %>%
  unnest(model_fit) %>%
  ungroup() %>%
  select(cancer, r_squared, adj_r_squared, p_value) %>%
  rename(model_p_value = p_value)
all_correlations <- all_expect_frequencies %>%
  group_by(cancer) %>%
  summarise(model_fit = list(
    calc_obs_pred_correlation(
      observed_allele_frequency,
      expected_allele_frequency
    )
  )) %>%
  unnest(model_fit)


#### ---- Statistics: Chi-Squared ---- ####


cancer_chisquared_res <- calc_chisquared_test(
  allele_df = alleles_for_each_cancer_obs_vs_pred,
  expected_freq_df = cancer_expect_frequencies
)
cancer_chisquared_res %>%
  filter(adj_p_value > 0.05) %>%
  select(cancer, kras_allele, p_value, adj_p_value) %>%
  knitr::kable(digits = 2)
# > |cancer |kras_allele | p_value| adj_p_value|
# > |:------|:-----------|-------:|-----------:|
# > |COAD   |G12A        |    0.69|        0.69|
# > |LUAD   |G12A        |    0.04|        0.06|
# > |LUAD   |G12D        |    0.16|        0.16|
# > |LUAD   |G12V        |    0.09|        0.12|
# > |MM     |G12A        |    0.58|        0.70|
# > |MM     |G12D        |    0.15|        0.41|
# > |MM     |G12R        |    0.61|        0.70|
# > |MM     |G12V        |    0.34|        0.55|
# > |MM     |G13D        |    0.05|        0.19|
# > |MM     |Q61L        |    1.00|        1.00|
# > |MM     |Q61R        |    0.29|        0.55|


all_chisquared_res <- calc_chisquared_test(
  allele_df = all_observed_alleles_obs_vs_pred,
  expected_freq_df = all_expect_frequencies
)
all_chisquared_res %>%
  filter(adj_p_value > 0.05) %>%
  select(cancer, kras_allele, p_value, adj_p_value) %>%
  knitr::kable(digits = 2)
# > |cancer |kras_allele | p_value| adj_p_value|
# > |:------|:-----------|-------:|-----------:|
# > |COAD   |G12C        |    0.11|        0.15|
# > |COAD   |G12R        |    0.14|        0.17|
# > |COAD   |G13D        |    0.50|        0.50|
# > |COAD   |Q61L        |    0.43|        0.46|
# > |MM     |G12A        |    0.53|        0.78|
# > |MM     |G12C        |    0.86|        0.92|
# > |MM     |G12D        |    0.56|        0.78|
# > |MM     |G12R        |    0.12|        0.24|
# > |MM     |G12V        |    0.61|        0.78|
# > |MM     |G13C        |    0.05|        0.11|
# > |MM     |G13D        |    0.76|        0.88|
# > |MM     |Q61L        |    0.38|        0.67|
# > |MM     |Q61R        |    0.92|        0.92|



#### ---- Plot: Distribution of probabilities ---- ####

return_one_followed_by_NA <- function(x) {
  c(x[[1]], rep(NA, length(x) - 1))
}


per_sample_allele_probability_plot <- all_kras_allele_predictions %>%
  group_by(cancer, kras_allele) %>%
  mutate(real_allele_freq = return_one_followed_by_NA(real_allele_freq)) %>%
  ungroup() %>%
  ggplot(aes(x = allele_prob)) +
  facet_grid(cancer ~ kras_allele, scales = "free_y") +
  geom_density(aes(color = kras_allele), fill = NA) +
  geom_vline(aes(xintercept = real_allele_freq, color = kras_allele),
    lty = 2, alpha = 0.6
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0.25, 0.75, length.out = 3),
    labels = scales::percent_format()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  scale_color_manual(values = short_allele_pal) +
  theme_bw(base_size = 7, base_family = "arial") +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 4.5, vjust = 1),
    axis.text.y = element_text(size = 4.5, hjust = 1),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "probability per sample",
    y = "density",
    title = "Per sample probability of each KRAS allele",
    subtitle = "Separated by cancer; color relates to KRAS allele; vertical lines are the actual frequency"
  )
ggsave_wrapper(
  per_sample_allele_probability_plot,
  plot_path(GRAPHS_DIR, "per_sample_allele_probability_plot.svg"),
  "wide"
)


#### ---- Plot: Predicted vs. Observed ---- ####

plot_kras_allele_predictions <- function(cancer,
                                         data,
                                         p_val_cut = 0.05,
                                         zero_axis_lines = FALSE) {
  pval_labels <- c(
    glue("p < {p_val_cut}"),
    glue("p ≥ {p_val_cut}")
  )

  codon_pal <- codon_palette[names(codon_palette) != "Other"]

  mod_data <- data %>%
    mutate(
      lower_ci = scales::squish(lower_ci, range = c(0, 1)),
      upper_ci = scales::squish(upper_ci, range = c(0, 1)),
      is_significant = ifelse(
        p_value < p_val_cut, pval_labels[[1]], pval_labels[[2]]
      ),
      is_significant = factor(is_significant, levels = pval_labels),
      codon = str_extract(kras_allele, "[:digit:]+|WT"),
      codon = factor(codon, levels = names(codon_pal))
    )

  max_val <- max(c(
    mod_data$observed_allele_frequency,
    mod_data$upper_ci
  ))

  p <- mod_data %>%
    ggplot(aes(
      x = expected_allele_frequency,
      y = observed_allele_frequency
    ))

  if (zero_axis_lines) {
    p <- p +
      geom_vline(xintercept = 0, size = 0.2, color = "grey50") +
      geom_hline(yintercept = 0, size = 0.2, color = "grey50")
  }

  p <- p +
    geom_abline(lty = 2, size = 0.6, color = "grey60") +
    geom_linerange(
      aes(
        xmin = lower_ci,
        xmax = upper_ci
      ),
      color = "grey30",
      alpha = 0.4,
      size = 0.4
    ) +
    geom_text_repel(
      aes(label = kras_allele),
      size = 2.2,
      family = "Arial",
      seed = 0,
      force = 1,
      box.padding = unit(2, "mm"),
      segment.alpha = 0.5,
      segment.color = "grey30",
      segment.size = 0.3,
      min.segment.length = unit(5, "mm"),
      max.time = 1,
      max.iter = 100000,
      max.overlaps = Inf
    ) +
    geom_point(
      aes(shape = is_significant, color = codon),
      size = 1.3
    ) +
    scale_color_manual(
      values = codon_pal,
      drop = FALSE,
      guide = guide_legend(order = 10)
    ) +
    scale_x_continuous(
      limits = c(0, max_val),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_y_continuous(
      limits = c(0, max_val),
      expand = expansion(mult = c(0, 0.03)),
      labels = function(x) {
        ifelse(x == 0, "", x)
      }
    ) +
    scale_shape_manual(
      values = c(17, 16),
      drop = FALSE,
      guide = guide_legend(
        title.position = "top",
        label.position = "top",
        order = 20
      )
    ) +
    coord_fixed() +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_markdown(hjust = 0.05),
      strip.background = element_blank(),
      legend.title = element_markdown(hjust = 0.5, vjust = 0.5)
    ) +
    labs(
      x = "predicted",
      y = "observed",
      shape = "χ<sup>2</sup> test",
      color = "codon",
      title = cancer
    )

  if (is.null(mod_data$cor_estimate)) {
    return(p)
  }

  cor_coef <- round(mod_data$cor_estimate[[1]], 3)
  cor_coef <- str_pad(cor_coef, width = 5, side = "right", pad = "0")
  subtitle <- glue("Pearson correlation: {cor_coef}")

  p <- p +
    labs(subtitle = subtitle)
  return(p)
}


save_kras_allele_predictions <- function(cancer, plt, gl_template,
                                         save_rds = TRUE) {
  ggsave_wrapper(
    plt,
    plot_path(GRAPHS_DIR, as.character(glue(gl_template))),
    size = "small"
  )

  if (save_rds) {
    saveFigRds(plt, as.character(glue(gl_template)))
  }

  invisible(plt)
}


cancer_expect_frequencies %>%
  left_join(cancer_chisquared_res, by = c("cancer", "kras_allele")) %>%
  left_join(cancer_rsquareds, by = c("cancer")) %>%
  left_join(cancer_correlations, by = c("cancer")) %>%
  group_by(cancer) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    plt = map2(cancer, data, plot_kras_allele_predictions),
    plt = map2(cancer, plt, save_kras_allele_predictions,
      gl_template = "{cancer}_predict-allele-freq_scatter.svg"
    )
  )

all_expect_frequencies %>%
  left_join(all_chisquared_res, by = c("cancer", "kras_allele")) %>%
  left_join(all_rsquareds, by = c("cancer")) %>%
  left_join(all_correlations, by = c("cancer")) %>%
  group_by(cancer) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    plt = map2(cancer, data, plot_kras_allele_predictions,
      zero_axis_lines = TRUE
    ),
    plt = map2(cancer, plt, save_kras_allele_predictions,
      gl_template = "{cancer}_predict-ALL-allele-freq_scatter.svg"
    )
  )


#### ---- Correlation coefficients of G12 mutations ---- ####

cancer_expect_frequencies %>%
  filter(str_detect(kras_allele, "G12")) %>%
  group_by(cancer) %>%
  summarise(model_fit = list(
    calc_obs_pred_correlation(
      observed_allele_frequency,
      expected_allele_frequency
    )
  )) %>%
  unnest(model_fit) %>%
  knitr::kable(digits = 3)
# > |cancer | cor_estimate| cor_p_value| cor_conf_low| cor_conf_high|cor_method                           |
# > |:------|------------:|-----------:|------------:|-------------:|:------------------------------------|
# > |COAD   |        0.759|       0.137|       -0.374|         0.983|Pearson's product-moment correlation |
# > |LUAD   |        0.870|       0.130|       -0.556|         0.997|Pearson's product-moment correlation |
# > |MM     |        0.972|       0.028|        0.167|         0.999|Pearson's product-moment correlation |
# > |PAAD   |        0.677|       0.323|       -0.813|         0.992|Pearson's product-moment correlation |
