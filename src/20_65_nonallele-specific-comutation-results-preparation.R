# Prepare results from non-allele-specific comutation analysis.

cache(
  "rc_test_nonallele_results",
  depends = "cancer_coding_muts_df",
  {
    rc_test_nonallele_results <- list.files(
      file.path("data", "rc-test-nonallelespec", "output", "exclusivity"),
      full.names = TRUE,
      pattern = "rds$"
    ) %>%
      unlist() %>%
      initial_rc_results_processing() %>%
      join_relevant_counts(
        muts_df = cancer_coding_muts_df
      )

    return(rc_test_nonallele_results)
  }
)
