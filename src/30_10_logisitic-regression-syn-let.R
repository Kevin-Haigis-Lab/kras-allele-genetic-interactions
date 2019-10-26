
#### ---- Multinomial logistic regression ---- ####

# I used multinomial logistic regression identify allele-specific synthetic
# lethality.


library(nnet)

#### ---- Multinomial logisitic regresssion ---- ####
# Regress on all of the alleles at once for each gene individually.

# Calculate the p-values using the z-test (Wald's).
# `m`: a model fitted using `nnet::multinom()`
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


# A wrapper around multinomial logisitic regression.
# It uses `nnet::multinom()` for fitting the model.
multinom_logisitc_regression_wrapper <- function(df) {
    df$allele <- factor(df$allele, levels = get_allele_factor_levels(df$allele))
    fit <- suppressMessages(nnet::multinom(allele ~ gene_effect, data = df))
    fit_pvals <- multinomial_logreg_pvalue(fit)
    return(list(
        logreg_model = fit,
        logreg_pvals = fit_pvals
    ))
}
multinom_logisitc_regression_wrapper <- memoise::memoise(multinom_logisitc_regression_wrapper)


# Conduct the multinomial modeling of the DepMap data.
log_model1_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        allele_logreg = purrr::map(data, multinom_logisitc_regression_wrapper)
    )




#### ---- Logisitic gression one each allele vs the others ---- ####
# Instead of regressing on all of the alleles at once, this modeling attempt
# regresses on each allele vs all other samples for each gene.
#
# I will use a Bonferroni correction within each gene (ie. divide the p-values
# for each gene by the number of tests).


# Prepare the data frame for comparing one allele `for_allele` against the rest.
prepare_log_model_data <- function(df, for_allele) {
    df %>%
        mutate(allele = ifelse(allele == !!for_allele, 1, 0))
}


# A wrapper for running logistic regression on all of the alleles, one at a time.
logisitc_regression_wrapper <- function(df) {
    # Fit a model for a single allele.
    log_wrap <- function(regression_allele, data) {
        mod_data <- prepare_log_model_data(data, regression_allele)
        glm(allele ~ gene_effect, family = "binomial", data = mod_data)
    }

    alleles <- sort(unique(df$allele))
    results <- tibble(allele = alleles[alleles != "WT"]) %>%
        mutate(fit_logreg = purrr::map(allele, log_wrap, data = df))
    return(results)
}


# Get the p-values for multiple fitted glm.
# With correct as `TRUE`, a Bonferroni correction is applied.
logreg_pvalues <- function(fit_glms, with_correction = TRUE) {
    pvals <- purrr::map_dbl(fit_glms, ~ broom::tidy(.x)$p.value[[2]])
    if (with_correction) pvals <- pvals * length(fit_glms)
    return(pvals)
}


# A plot of the logistic regression between one allele and the rest.
logreg_plot <- function(cancer, hugo_symbol, allele, ...) {
    p <- model_data %>%
        filter(cancer == !!cancer & hugo_symbol %in% !!hugo_symbol) %>%
        mutate(y = ifelse(allele == !!allele, 1, 0)) %>%
        ggplot(aes(x = gene_effect, y = y)) +
        geom_point(aes(color = allele)) +
        scale_color_manual(values = short_allele_pal) +
        geom_smooth(method = "glm", method.args = list(family = "binomial")) +
        theme_bw(base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5)
        ) +
        labs(
            x = "dependency",
            y = glue("allele (other: 0, {allele}: 1)"),
            title = hugo_symbol
        )
    return(p)
}


# Make and save a logistic regression plots.
save_logreg_plot <- function(cancer, hugo_symbol, allele, ...) {
    p <- logreg_plot(cancer, hugo_symbol, allele)
    save_path <- plot_path("30_10_logisitic-regression-syn-let_model2",
                           glue("{cancer}_{allele}_{hugo_symbol}_logreg.svg"))
    ggsave_wrapper(p, save_path, "medium")
}


# Run logistic regression on the KRAS alleles.
log_model2_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        allele_logreg = purrr::map(data, logisitc_regression_wrapper)
    ) %>%
    select(-data) %>%
    unnest(allele_logreg) %>%
    group_by(hugo_symbol, cancer) %>%
    mutate(
        p_val = logreg_pvalues(fit_logreg, with_correction = FALSE),
        p_val_corrected = p_val * length(p_val)
    ) %>%
    ungroup()


# Show number of hits for each cancer.
log_model2_tib %>%
    filter(p_val_corrected < 0.05) %>%
    count(cancer)



# Save the few for COAD and PAAD
log_model2_tib %>%
    filter(p_val_corrected < 0.05 & cancer != "LUAD") %>%
    pwalk(save_logreg_plot)


# Save the top 10 for LUAD
luad_plots <- log_model2_tib %>%
    filter(p_val_corrected < 0.05 & cancer == "LUAD") %>%
    top_n(9, -p_val) %>%
    select(cancer, hugo_symbol, allele) %>%
    pmap(logreg_plot)
luad_plot <- cowplot::plot_grid(plotlist = luad_plots, nrow = 3)
save_path <- plot_path("30_10_logisitic-regression-syn-let_model2",
                       "LUAD_top9_logreg.svg")
cowplot::save_plot(save_path, luad_plot, base_height = 12, base_width = 12)

