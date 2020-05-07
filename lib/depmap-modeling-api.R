# An API for common manipulations of `depmap_model_workflow_res`: the results
# of the DepMap modeling workflow.

#### ---- Filtering the workflow ---- ####

# A standard way to unnest the RNA-expr. linear model results.
unnest_rna_lm_depmap_model_workflow <- function(df) {
    df %>%
        unnest(rna_lm_stats,  names_sep = "_")
}

# A single function that will filter the results of the above workflow.
# This will maintain consistency throughout the project.
filter_depmap_model_workflow_res <- function(df) {
    df %>%
        unnest_rna_lm_depmap_model_workflow() %>%
        filter(rna_lm_stats_p_value >= 0.05) %>%
        filter(ismut_pvalue >= 0.05) %>%
        filter(allele_aov_pvalue < 0.01) %>%
        filter(ova_any_sig)
}
