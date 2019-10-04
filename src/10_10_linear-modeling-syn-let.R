

#### ---- Linear Model 1 ---- ####
info(logger, "Beginning linear model 1.")

get_allele_factor_levels <- function(als) {
    if (any(als == "WT")) {
        all_but_wt <- als[als != "WT"] %>% sort()
        return(c("WT", all_but_wt))
    } else {
        return(sort(als))
    }
}


lm1_prepare_data <- function(tib) {
    allele_levels <- get_allele_factor_levels(unique(tib$allele))
    mod_tib <- tib %>%
        mutate(
            allele = factor(allele, levels = allele_levels),
            is_altered = as.numeric(is_altered),
            rna_scaled = scale(rna_expression)[, 1]
        ) %>%
        select(gene_effect, allele, is_altered, rna_scaled)
    if (all(is.na(mod_tib$rna_scaled))) { mod_tib$rna_scaled <- 0 }
    return(mod_tib)
}


lm1_function <- function(hugo_symbol, cancer, data, ...) {
    data_mod <- lm1_prepare_data(data)
    fit <- try(lm(
        gene_effect ~ allele + is_altered + rna_scaled,
        data = data_mod
    ))
    fit_aic <- MASS::stepAIC(fit, direction = "both", trace = FALSE)
    return(fit_aic)
}


test_lm <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    slice(1:100) %>% # TEST
    pmap(lm1_function)
