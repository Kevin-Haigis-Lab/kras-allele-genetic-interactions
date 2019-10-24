
#### ---- Multinomial logistic regression ---- ####

# I used multinomial logistic regression identify allele-specific synthetic
# lethality.


library(nnet)

multinomial_logreg_pvalue <- function(m) {
    zvals <- summary(m)$coefficients / summary(m)$standard.errors
    pvals <- (1 - pnorm(abs(zvals), 0, 1)) * 2
    return(pvals)
}



# put the alleles in the correct order of turning into factors
#   if there is WT, put that first
#   all of the other alleles are in alphabetical order
get_allele_factor_levels <- function(als) {
    als <- unique(als)
    if (any(als == "WT")) {
        all_but_wt <- sort(als[als != "WT"])
        return(c("WT", all_but_wt))
    } else {
        return(sort(als))
    }
}




logisitc_regression_wrapper <- function(df) {
    df$allele <- factor(df$allele, levels = get_allele_factor_levels(df$allele))
    fit <- suppressMessages(nnet::multinom(allele ~ gene_effect, data = df))
    fit_pvals <- multinomial_logreg_pvalue(fit)
    return(list(
        logreg_model = fit,
        logreg_pvals = fit_pvals
    ))
}



# conduct the first modeling attempt
model2_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    ungroup() %>%
    sample_n(1e3) %>%
    mutate(
        allele_logreg = purrr::map(data, logisitc_regression_wrapper)
    )

# info(logger, "Caching results of model 1.")
# cache("model1_tib")




# TODO:
#   - logistic regression on each allele vs. others