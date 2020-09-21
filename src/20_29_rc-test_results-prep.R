# Prepare results from RC test.
cache(
  "rc_test_results",
  depends = "cancer_coding_muts_df",
  {
    rc_test_results <- purrr::map(
      c("comutation", "exclusivity"),
      ~ list(list.files(
        file.path("data", "rc-test", "output", .x),
        full.names = TRUE,
        pattern = "rds"
      ))
    ) %>%
      unlist() %>%
      initial_rc_results_processing() %>%
      join_relevant_counts(
        muts_df = cancer_coding_muts_df
      )

    return(rc_test_results)
  }
)
