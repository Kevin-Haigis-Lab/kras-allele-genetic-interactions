# Analyze the results of the non-allele-specific comutation analysis.



source(file.path("src", "20_34_rc-fisher-assessment-processes.R"))

rc_test_nonallele_results %>%
  select(-rc_test_type, -allele) %>%
  rename(rc_test_type = kras_allele) %>%
  select(
    cancer, hugo_symbol, p_val, t_AM, t_BM_ge,
    num_samples_per_cancer:num_mut_per_cancer_allele
  ) %>%
  assess_rc_test_significance()


nonallele_specific_increased_comutation_df %>%
  #> do some prep here %>%
  assess_fisher_test_significance()