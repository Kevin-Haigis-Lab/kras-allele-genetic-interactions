# The hits from the survival analysis that are worth following up on.

survival_analysis_hits <- tibble::tribble(
    ~ cancer, ~ allele, ~ hugo_symbol, ~ comutation_interaction,
      "LUAD",   "G12C",      "ARID1A",                "reduced",
      "LUAD",   "G12C",      "CHRNB4",              "increased",
      "LUAD",   "G12C",       "VN1R2",              "increased",
      "LUAD",   "G12C",      "ZNF445",              "increased",
      "LUAD",   "G12C",     "ZNF804A",                "reduced"
)


# Check for correct labels in `comutation_interaction` column.
CHECK1 <- survival_analysis_hits %>%
    filter(!comutation_interaction %in% c("reduced", "increased"))
if (nrow(CHECK1) > 0) {
    stop("Incorrect labels for comutation interaction.")
}
