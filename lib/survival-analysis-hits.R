# The hits from the survival analysis that are worth following up on.

survival_analysis_hits <- tibble::tribble(
     ~cancer,  ~interaction_allele,  ~hugo_symbol,  ~comutation_interaction, ~likelihood_ratio_test_pval, ~allele_pval, ~comutation_pval,
      "LUAD",               "G12C",      "ARID1A",                "reduced",                          NA,           NA,               NA,
      "LUAD",               "G12C",      "CHRNB4",              "increased",                        0.04,         0.17,           0.0288,
      "LUAD",               "G12C",       "VN1R2",              "increased",                        0.02,       0.0661,           0.0161,
      "LUAD",               "G12C",      "ZNF445",              "increased",                        0.03,       0.0542,           0.0244,
      "LUAD",               "G12C",     "ZNF804A",                "reduced",                        0.04,       0.0705,           0.0848,
)


# Check for correct labels in `comutation_interaction` column.
CHECK1 <- survival_analysis_hits %>%
    filter(!comutation_interaction %in% c("reduced", "increased"))
if (nrow(CHECK1) > 0) {
    stop("Incorrect labels for comutation interaction.")
}
