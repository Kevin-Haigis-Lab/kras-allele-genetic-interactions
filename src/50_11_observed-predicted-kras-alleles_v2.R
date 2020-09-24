# Observed vs. predicted KRAS allele frequency.
# This is another analysis with the same goal. I just want to make sure
#   it was done correctly. Hopefully, I get the same results as the first time.

GRAPHS_DIR <- "50_11_observed-predicted-kras-alleles_v2"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)

mut_sig_descriptions <- mutational_signatures_df %>%
  distinct(signature, description)

artifact_signatures <- mut_sig_descriptions %>%
  filter(description == "artifact") %>%
  mutate(signature = paste0("sig", signature)) %>%
  pull(signature)

msg_sigs <- paste0(artifact_signatures, collapse = ", ")
message(glue("The following are artifact signatures: {msg_sigs}"))


oncogenic_alleles <- unique(kras_trinucleotide_contexts$kras_allele)
message("The following are the oncogenic KRAS alleles:")
print(oncogenic_alleles)


#### ---- KRAS allele frequency ---- ####

# Calculate the frequency of the KRAS mutant alleles.
# This only includes the oncogenic KRAS alleles.
calc_frequency_of_oncogenic_kras <- function(alleles) {
  allele_tbl <- table(alleles)
  allele_tbl <- allele_tbl[names(allele_tbl) %in% oncogenic_alleles]
  allele_tbl <- as_tibble(allele_tbl) %>%
    mutate(n = n / sum(n)) %>%
    set_names(c("kras_allele", "frequency"))

  missing_idx <- !(oncogenic_alleles %in% allele_tbl$kras_allele)
  missing_alleles <- oncogenic_alleles[missing_idx]

  if (length(missing_alleles) > 0) {
    missing_tbl <- tibble(
      kras_allele = missing_alleles,
      frequency = 0.0
    )
    allele_tbl <- bind_rows(allele_tbl, missing_tbl)
  }
  return(allele_tbl)
}


real_kras_allele_freq <- cancer_full_coding_muts_df %>%
  filter(cancer != "SKCM") %>%
  distinct(cancer, tumor_sample_barcode, ras_allele) %>%
  mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
  group_by(cancer) %>%
  summarise(
    kras_freq = list(calc_frequency_of_oncogenic_kras(kras_allele))
  ) %>%
  unnest(kras_freq) %>%
  dplyr::rename(real_allele_freq = frequency)



#### ---- Matrix A: KRAS allele x mutational signature ---- ####

kras_by_mutsig_matrix <- kras_trinucleotide_contexts %>%
  left_join(mutational_signature_spectra, by = "tricontext") %>%
  select(kras_allele, signature, composition) %>%
  mutate(signature = paste0("sig", signature)) %>%
  filter(!signature %in% artifact_signatures) %>%
  group_by(kras_allele, signature) %>%
  summarise(composition = sum(composition)) %>%
  ungroup() %>%
  pivot_wider(kras_allele,
    names_from = signature,
    values_from = composition
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  column_to_rownames("kras_allele") %>%
  as.matrix()

kras_by_mutsig_matrix[, 1:5]

# Add a row to take up the remainder of the mutational signature.
residual_signature <- matrix(1 - colSums(kras_by_mutsig_matrix), nrow = 1)
rownames(residual_signature) <- "remainder"
kras_by_mutsig_matrix <- rbind(kras_by_mutsig_matrix, residual_signature)

kras_by_mutsig_matrix[, 1:5]

# (Plots of these values per allele are at the end.)


#### ---- Matrix B: mutational signature x tumor sample barcode ---- ####

MIN_NUM_MUTATION_PER_SAMPLE <- 10

REMOVE_TSB <- trinucleotide_mutations_df %>%
  count(tumor_sample_barcode) %>%
  filter(n < MIN_NUM_MUTATION_PER_SAMPLE) %>%
  pull(tumor_sample_barcode)

message(glue(
  "Removing {length(REMOVE_TSB)} samples because they have too few mutations"
))


mutsig_by_tsb_matrix <- mutational_signatures_df %>%
  filter(!tumor_sample_barcode %in% REMOVE_TSB) %>%
  select(tumor_sample_barcode, signature, contribution) %>%
  mutate(signature = paste0("sig", signature)) %>%
  filter(!signature %in% artifact_signatures) %>%
  pivot_wider(tumor_sample_barcode,
    names_from = signature,
    values_from = contribution
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  column_to_rownames("tumor_sample_barcode") %>%
  as.matrix() %>%
  t()

mutsig_by_tsb_matrix[, 1:5]

# Remove samples with no mutational signature information.
# (All of their mutations were contributed to artifact).
sum(colSums(mutsig_by_tsb_matrix) == 0)
mutsig_by_tsb_matrix <- mutsig_by_tsb_matrix[, colSums(mutsig_by_tsb_matrix) > 0]
all(colSums(mutsig_by_tsb_matrix) > 0)


#### ---- Organize matrices A and B ---- ####

# Set column names of A equal to row names of B.
all_sigs <- unique(paste0("sig", mutational_signature_spectra$signature))
all_sigs <- setdiff(all_sigs, artifact_signatures)
all(all_sigs %in% colnames(kras_by_mutsig_matrix))
all(all_sigs %in% rownames(mutsig_by_tsb_matrix))

kras_by_mutsig_matrix <- kras_by_mutsig_matrix[, all_sigs]
mutsig_by_tsb_matrix <- mutsig_by_tsb_matrix[all_sigs, ]

identical(colnames(kras_by_mutsig_matrix), rownames(mutsig_by_tsb_matrix))

kras_by_mutsig_matrix[1:5, 1:5]
mutsig_by_tsb_matrix[1:5, 1:5]


#### ---- Matrix C: KRAS allele x tumor sample barcode ---- ####

kras_by_tsb_matrix <- kras_by_mutsig_matrix %*% mutsig_by_tsb_matrix

kras_by_tsb_matrix[, 1:5]
colSums(kras_by_tsb_matrix[, 1:5])


#### ---- Analysis of results ---- ####

# Tumor sample information to add back to results.
sample_information <- trinucleotide_mutations_df %>%
  distinct(cancer, tumor_sample_barcode, target, kras_allele) %>%
  dplyr::rename(real_kras_allele = kras_allele)


# Matrix to long tibble.
kras_by_tsb_df <- kras_by_tsb_matrix %>%
  as.data.frame() %>%
  rownames_to_column("kras_allele") %>%
  as_tibble() %>%
  pivot_longer(-kras_allele,
    names_to = "tumor_sample_barcode",
    values_to = "probability"
  ) %>%
  left_join(sample_information, by = "tumor_sample_barcode")


#### ---- Allele frequency as most likely per sample ---- ####
# To predict the frequency of the KRAS alleles, we can assign each sample its
# most likely KRAS allele and then find the frequency of these
# predicted alleles.

# without remainder
kras_freq_per_sample <- kras_by_tsb_df %>%
  filter(kras_allele != "remainder") %>%
  group_by(tumor_sample_barcode, cancer) %>%
  summarise(pred_kras_allele = kras_allele[which.max(probability)]) %>%
  ungroup() %>%
  count(cancer, pred_kras_allele, name = "pred_count") %>%
  right_join(real_kras_allele_freq,
    by = c("pred_kras_allele" = "kras_allele", "cancer")
  ) %>%
  mutate(pred_count = ifelse(is.na(pred_count), 0, pred_count)) %>%
  group_by(cancer) %>%
  mutate(pred_allele_freq = pred_count / sum(pred_count))


kras_freq_per_sample_scatter <- kras_freq_per_sample %>%
  ggplot(aes(x = real_allele_freq, y = pred_allele_freq)) +
  facet_wrap(. ~ cancer, scales = "free") +
  geom_abline(lty = 2, color = "grey80") +
  geom_point(size = 0.7, color = "grey30") +
  ggrepel::geom_text_repel(
    aes(label = pred_kras_allele),
    size = 2, family = "Arial",
    segment.size = 0.1, segment.alpha = 0.5, segment.color = "grey50",
    seed = 0
  ) +
  theme_bw(base_size = 7, base_family = "arial") +
  theme(
    strip.background = element_blank()
  ) +
  labs(x = "observed", y = "predicted")
ggsave_wrapper(
  kras_freq_per_sample_scatter,
  plot_path(GRAPHS_DIR, "kras_freq_per_sample_scatter.svg"),
  "medium"
)


#### ---- Allele frequency averaged across samples ---- ####

# The minimum number of samples with a KRAS allele for the KRAS allele to be
# included in the calculation (per cancer).

kras_freq_sample_avg <- kras_by_tsb_df %>%
  filter(kras_allele != "remainder") %>%
  group_by(tumor_sample_barcode) %>%
  mutate(probability = probability / sum(probability)) %>%
  group_by(cancer, kras_allele) %>%
  summarise(
    avg_probability = mean(probability),
    sd_probability = sd(probability)
  ) %>%
  ungroup() %>%
  left_join(real_kras_allele_freq,
    by = c("cancer", "kras_allele")
  )


kras_freq_sample_avg %>%
  group_by(cancer) %>%
  summarise(total_probability = sum(avg_probability))


plot_obs_pred_pointrange_scatter <- function(cancer, data) {
  max_val <- max(c(data$ymax, data$real_allele_freq))
  data$ymin <- map_dbl(data$ymin, ~ max(.x, 0))

  ggplot(data, aes(x = real_allele_freq, y = avg_probability)) +
    geom_abline(lty = 2, color = "grey80") +
    geom_pointrange(
      aes(ymin = ymin, ymax = ymax),
      fatten = 1.6, size = 0.4, color = "grey30"
    ) +
    ggrepel::geom_text_repel(
      aes(label = kras_allele),
      size = 4, family = "Arial",
      segment.size = 0.7, segment.alpha = 0.5, segment.color = "grey50",
      point.padding = unit(2, "mm"),
      min.segment.length = unit(5, "mm"),
      seed = 0
    ) +
    theme_bw(base_size = 10, base_family = "arial") +
    theme(
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(
      limits = c(0, max_val),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, max_val),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_fixed() +
    labs(x = "observed", y = "predicted", title = cancer)
}


# Save the plot.
save_kras_freq_pointrange_scatter <- function(cancer, plt, gl_template) {
  ggsave_wrapper(
    plt,
    plot_path(GRAPHS_DIR, glue(gl_template)),
    "small"
  )
  invisible(plt)
}


# Write the data to file.
write_kras_frequencies <- function(cancer, data, gl_template) {
  write_tsv(
    data,
    table_path(GRAPHS_DIR, glue(gl_template))
  )
  invisible(data)
}

# All oncogenic KRAS alleles for every cancer.
kras_freq_sample_avg %>%
  mutate(
    ymin = avg_probability - sd_probability,
    ymax = avg_probability + sd_probability
  ) %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(
    plt = map2(cancer, data, plot_obs_pred_pointrange_scatter),
    plt = map2(
      cancer, plt, save_kras_freq_pointrange_scatter,
      gl_template = "{cancer}_kras-freq_sample-avg_all-alleles.svg"
    ),
    data = map2(
      cancer, data, write_kras_frequencies,
      gl_template = "{cancer}_kras-freq_sample-avg_all-alleles.tsv"
    )
  )


# Only the KRAS alleles above a certain frequency in each cancer.
# The cut-off is a 2% frequency of the KRAS allele in the cancer.
kras_freq_sample_avg %>%
  filter(real_allele_freq >= 0.02) %>%
  group_by(cancer) %>%
  mutate(
    real_allele_freq = real_allele_freq / sum(real_allele_freq),
    sd_probability = sd_probability / sum(avg_probability),
    avg_probability = avg_probability / sum(avg_probability)
  ) %>%
  ungroup() %>%
  mutate(
    ymin = avg_probability - sd_probability,
    ymax = avg_probability + sd_probability
  ) %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(
    plt = map2(cancer, data, plot_obs_pred_pointrange_scatter),
    plt = map2(
      cancer, plt, save_kras_freq_pointrange_scatter,
      gl_template = "{cancer}_kras-freq_sample-avg_cancer-alleles.svg"
    ),
    data = map2(
      cancer, data, write_kras_frequencies,
      gl_template = "{cancer}_kras-freq_sample-avg_cancer-alleles.tsv"
    )
  )


#### ---- Plot the contribution of the mut sigs to KRAS alleles ---- ####

mutsig_contribution_to_alleles <- kras_by_mutsig_matrix %>%
  as.data.frame() %>%
  rownames_to_column("kras_allele") %>%
  as_tibble() %>%
  pivot_longer(-kras_allele,
    names_to = "signature",
    values_to = "contribution"
  ) %>%
  mutate(signature = str_remove(signature, "^sig")) %>%
  left_join(mut_sig_descriptions, by = "signature") %>%
  group_by(kras_allele, description) %>%
  summarise(contribution = sum(contribution)) %>%
  ungroup() %>%
  mutate(description = factor(description,
    levels = names(mutsig_descrpt_pal)
  )) %>%
  ggplot(aes(x = kras_allele, y = contribution, fill = description)) +
  geom_col(position = "fill") +
  scale_fill_manual(
    values = mutsig_descrpt_pal,
    guide = guide_legend(nrow = 2)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_markdown(hjust = 0.5)
  ) +
  labs(
    fill = "mut. sig.",
    title = "Contribution of the mutational signatures to the *KRAS* alleles"
  )
ggsave_wrapper(
  mutsig_contribution_to_alleles,
  plot_path(GRAPHS_DIR, "mutsig_contribution_to_alleles.svg"),
  "wide"
)
