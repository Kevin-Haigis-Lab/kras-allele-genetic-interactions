# Methods used in predicting KRAS allele frequencies from mutational signatures.

#### ---- Utilities ---- ####


check_cached_data_exists <- function(x) {
  if (!exists(x)) {
    stop(glue::glue("data frame missing: `{x}`"))
  }
  invisible(TRUE)
}

#### ---- Data getters ---- ####


# Get a data frame of mutational signature ID and descriptive name.
get_mut_sig_descriptions <- function() {
  check_cached_data_exists("mutational_signatures_df")
  mutational_signatures_df %>%
    distinct(signature, description)
}


# List of artifact signatures.
get_artifact_signatures <- function() {
  check_cached_data_exists("mutational_signatures_df")
  mut_sig_descriptions %>%
    filter(description == "artifact") %>%
    mutate(signature = paste0("sig", signature)) %>%
    pull(signature)
}


#### ---- Artifact signatures ---- ####


# SUm up all artifact signature contributions.
total_artifact_contribution <- function(mut_sig_df, artifact_sigs) {
  mut_sig_df %>%
    filter(signature %in% !!artifact_sigs) %>%
    group_by(tumor_sample_barcode) %>%
    summarise(contribution = sum(contribution)) %>%
    ungroup()
}

#### ---- Allele-frequency tables ---- ####


# Calculate the frequency of the KRAS mutant alleles.
# This only includes the oncogenic KRAS alleles.
calc_frequency_of_oncogenic_kras <- function(alleles, oncogenic_alleles) {
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


#### ---- Calculate predicted frequency ---- ####


count_tricontext_allele_mutations <- function(tricontext_mut_data,
                                              allele_df,
                                              tricontext_counts_df,
                                              kras_tricontexts,
                                              real_allele_freq,
                                              remove_samples = NULL,
                                              remove_cancers = c("SKCM")) {
  tricontext_mut_data %>%
    filter(!(cancer %in% !!remove_cancers)) %>%
    filter(!(tumor_sample_barcode %in% !!remove_samples)) %>%
    filter(hugo_symbol != "KRAS") %>%
    count(
      cancer, tumor_sample_barcode, target, context, tricontext,
      name = "tumor_count"
    ) %>%
    left_join(
      tricontext_counts_df,
      by = c("context", "target")
    ) %>%
    mutate(tumor_count_norm = tumor_count / genome_count) %>%
    inner_join(
      kras_tricontexts,
      by = c("context", "tricontext")
    ) %>%
    group_by(target) %>%
    complete(
      nesting(cancer, tumor_sample_barcode, target),
      nesting(kras_allele, kras_codon, context, tricontext, genome_count),
      fill = list(
        tumor_count = 0,
        tumor_count_norm = 0
      )
    ) %>%
    ungroup() %>%
    right_join(
      allele_df,
      by = c("cancer", "kras_allele")
    ) %>%
    group_by(cancer, tumor_sample_barcode, target, kras_allele) %>%
    summarise(allele_prob = sum(tumor_count_norm)) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    mutate(allele_prob = allele_prob / sum(allele_prob)) %>%
    ungroup() %>%
    left_join(
      real_allele_freq,
      by = c("cancer", "kras_allele")
    )
}

# Remove tumor samples with no mutations of the same "type" as the
# KRAS alleles. They effectively have no predictions.
filter_samples_with_no_predictions <- function(df) {
  df %>%
    group_by(cancer, tumor_sample_barcode) %>%
    filter(!all(is.na(allele_prob))) %>%
    ungroup()
}


#### ---- Statistics: boostrapping 95% CI ---- ####

# Calculate the expect frequency for the data of a single cancer.
calc_expected_frequency <- function(df) {
  df %>%
    group_by(kras_allele) %>%
    summarise(
      expected_allele_frequency = mean(allele_prob),
      expected_allele_frequency_lower25 = quantile(allele_prob, 0.25),
      expected_allele_frequency_lower75 = quantile(allele_prob, 0.75),
      observed_allele_frequency = unique(real_allele_freq)
    ) %>%
    ungroup()
}


# Bootstrap CIs for the expected frequencies.
calc_expected_frequency_boot <- function(data,
                                         index = NULL,
                                         all_alleles = NULL) {
  exp_freqs <- data[index, ] %>%
    unnest(data) %>%
    calc_expected_frequency() %>%
    select(-observed_allele_frequency)

  allele_results <- tibble(kras_allele = all_alleles) %>%
    left_join(exp_freqs, by = "kras_allele") %>%
    mutate(expected_allele_frequency = ifelse(
      is.na(expected_allele_frequency), 0, expected_allele_frequency
    )) %>%
    select(kras_allele, expected_allele_frequency)

  return(deframe(allele_results))
}


# Bootstrapping
boot_cancer_expect_frequncies <- function(cancer, df, R = 1e3) {
  nested_df <- df %>%
    group_by(tumor_sample_barcode) %>%
    nest()
  boot_res <- boot::boot(
    nested_df,
    calc_expected_frequency_boot,
    R = R,
    all_alleles = unique(df$kras_allele)
  )
  return(list(
    boot_obj = boot_res,
    alleles = unique(df$kras_allele)
  ))
}


# Extract the CIs from a boot object at an index.
get_conf_intervals <- function(boot_obj, conf, index) {
  boot::boot.ci(
    boot_obj,
    index = index,
    conf = conf,
    type = "perc"
  )$perc[, c(4, 5)]
}

# Extract 95% CI for each index from bootstrap results.
extract_boot_results <- function(boot_res, conf = 0.95) {
  f <- function(x, i) {
    ci <- get_conf_intervals(boot_res$boot_obj, conf = conf, index = i)
    tibble(kras_allele = x, lower_ci = ci[[1]], upper_ci = ci[[2]])
  }
  imap(boot_res$alleles, f) %>%
    bind_rows()
}


#### ---- Check calculations ---- ####

# Check that the sum of probabilities for each TSB is 1.
check_sum_of_probabilities <- function(d, prob_col, value = 1) {
  sum_probs <- d %>%
    summarise(sum_prob = sum({{ prob_col }})) %>%
    pull(sum_prob)
  return(all(near(sum_probs, value)))
}


#### ---- Statistics: R-squared ---- ####

calc_obs_pred_rsquared <- function(x, y) {
  dat <- tibble(x = x, y = y)
  lm(y ~ x, data = dat) %>%
    broom::glance() %>%
    janitor::clean_names()
}


calc_obs_pred_correlation <- function(x, y) {
  res <- cor.test(x, y) %>%
    broom::tidy() %>%
    janitor::clean_names() %>%
    select(estimate, p_value, conf_low, conf_high, method)
  colnames(res) <- paste0("cor_", colnames(res))
  return(res)
}


#### ---- Statistics: Chi-Squared ---- ####

# Chi-squared test to nest null hypothesis that observed and predicted
# frequency of an allele are the same.
allele_pred_obs_chisquared <- function(num_mut, num_tot, pred_freq) {
  pred_alleles <- num_tot * pred_freq
  mat <- matrix(
    c(
      num_mut, num_tot - num_mut,
      pred_alleles, num_tot - pred_alleles
    ),
    nrow = 2,
    byrow = TRUE
  )

  chisq.test(mat) %>%
    broom::glance() %>%
    janitor::clean_names()
}


calc_chisquared_test <- function(allele_df, expected_freq_df) {
  # Check for global data frames.
  if (!exists("cancer_full_coding_muts_df")) {
    print("data frame missing: `cancer_full_coding_muts_df`")
  }

  cancer_full_coding_muts_df %>%
    mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
    right_join(allele_df, by = c("cancer", "kras_allele")) %>%
    group_by(cancer, kras_allele) %>%
    summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer) %>%
    mutate(num_cancer_samples = sum(num_allele_samples)) %>%
    ungroup() %>%
    right_join(expected_freq_df, by = c("cancer", "kras_allele")) %>%
    group_by(cancer, kras_allele) %>%
    mutate(
      chi_squared_test = list(
        allele_pred_obs_chisquared(
          num_allele_samples,
          num_cancer_samples,
          expected_allele_frequency
        )
      )
    ) %>%
    ungroup() %>%
    select(cancer, kras_allele, chi_squared_test) %>%
    unnest(chi_squared_test) %>%
    group_by(cancer) %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
}
