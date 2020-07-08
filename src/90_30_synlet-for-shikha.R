# List of genes with synthetic lethal interactions with KRAS.
# comparisons:
#   1. syn. let. with KRAS (mut vs. WT)
#   2. syn. let. with G12D
#   3. syn. let. with G13D
#   4. syn. let. with G12V

library(effectsize)

GRAPHS_DIR <- "90_30_synlet-for-shikha"
reset_graph_directory(GRAPHS_DIR)

TABLES_DIR <- GRAPHS_DIR
reset_table_directory(TABLES_DIR)


set.seed(0)

depmap_modelling_df

coad_sl_res <- depmap_model_workflow_res %>%
    filter(cancer == "COAD") %>%
    filter_depmap_model_workflow_res()

coad_slcomut_res <- synlet_comut_model_res %>%
    filter(cancer == "COAD")


coad_dep_map_ids <- synlet_data %>%
    unnest(data) %>%
    u_pull(dep_map_id)

# Write out table of KRAS mutations in CCLE for Shikha.
ccle_cell_lines %>%
    filter(!is.na(cancer)) %>%
    left_join(ccle_kras_muts, by = "dep_map_id") %>%
    select(dep_map_id, ccle_name, stripped_cell_line_name,
           cancer, lineage, lineage_subtype, lineage_sub_subtype,
           lineage_molecular_subtype, disease, disease_subtype,
           primary_or_metastasis,
           sex, age, culture_medium,
           kras_allele = allele,
           kras_codon = codon,
           kras_cn = copy_number) %>%
    distinct() %>%
    mutate(kras_allele = ifelse(is.na(kras_allele), "WT", kras_allele),
           kras_codon = ifelse(is.na(kras_codon), "WT", kras_codon),
           kras_copynumber = ifelse(is.na(kras_cn), 2, kras_cn)) %>%
    filter(cancer == "COAD") %>%
    write_tsv(table_path(TABLES_DIR, "CCLE-cell-line-information.tsv"))



#### ---- Narrow down COAD KRAS allele-specific SL targets ---- ####

make_pretty_ci <- function(low, high) {
    low <- round(low, 2)
    high <- round(high, 2)
    glue("[{low}, {high}]")
}

calculate_allele_hedges_g <- function(allele, data) {
    mod_data <- data %>%
        mutate(is_allele = kras_allele == !!allele)
    hedges_g(gene_effect ~ is_allele, data = mod_data) %>%
        as_tibble() %>%
        janitor::clean_names() %>%
        mutate(hedges_g_ci = make_pretty_ci(ci_low, ci_high)) %>%
        select(hedges_g, hedges_g_ci) %>%
        mutate(hedges_g_interpret = interpret_d(hedges_g))
}

calculate_anova_omega_squared <- function(aov_mdl) {
    omega_squared(aov_mdl) %>%
        as_tibble() %>%
        janitor::clean_names() %>%
        mutate(omega_sq_ci = make_pretty_ci(ci_low, ci_high)) %>%
        select(omega_sq_partial, omega_sq_ci) %>%
        mutate(omega_sq_interpret = interpret_omega_squared(omega_sq_partial))
}


get_elastic_net_coefficients <- function(fit) {
    fit$elastic_model %>%
        coef() %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        as_tibble() %>%
        set_names(c("name", "coef_value")) %>%
        filter(abs(coef_value) > 0)
}


get_allele_estimated_effect <- function(allele, data) {
    mod_data <- data %>% mutate(is_allele = kras_allele == !!allele)
    lm(gene_effect ~ is_allele, data = mod_data) %>%
        broom::tidy() %>%
        select(term, estimate) %>%
        mutate(term = c("intercept_est", "allele_est")) %>%
        pivot_wider(names_from = term, values_from = estimate)
}



data_for_shiny_app <- coad_slcomut_res %>%
    left_join(coad_sl_res %>% select(cancer, hugo_symbol, allele_aov),
              by = c("cancer", "hugo_symbol")) %>%
    mutate(
        allele_hedges_g = map2(allele, data, calculate_allele_hedges_g),
        anova_omega_squared = map(allele_aov, calculate_anova_omega_squared),
        en_coefs = map(fit, get_elastic_net_coefficients),
        num_fit_coefs = map_dbl(en_coefs, ~ nrow(.x)),
        allele_in_final_model = map_lgl(en_coefs, ~ "kras_allele" %in% .x$name),
        allele_effect = map2(allele, data, get_allele_estimated_effect)
    ) %>%
    unnest(allele_hedges_g) %>%
    unnest(anova_omega_squared) %>%
    unnest(allele_effect) %>%
    mutate(intercept_allele_diff = abs(intercept_est - allele_est)) %>%
    select(hugo_symbol, allele, data, fit, num_coefs,
           allele_aov:intercept_allele_diff)

saveRDS(
    data_for_shiny_app,
    plot_path(GRAPHS_DIR, "data_for_shiny_app.rds")
)
