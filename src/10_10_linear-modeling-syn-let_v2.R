################################################################################
################################################################################

#### ---- Just for local work, not needed on O2 ---- ####

library(glue)
library(broom)
library(ggfortify)
library(patchwork)
library(jhcutils)
library(magrittr)
library(ggpubr)
library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")

filter_depmap_by_allele_count <- function(df, min = 3) {
    count_tib <- df %>%
        distinct(cancer, dep_map_id, kras_allele) %>%
        add_count(cancer, kras_allele, name = "n") %>%
        filter(n >= !!min) %>%
        select(cancer, dep_map_id, kras_allele)

    df %>%
        inner_join(count_tib, by = c("cancer", "dep_map_id", "kras_allele"))
}


load("cache/depmap_modelling_df.RData")
depmap_modelling_df %<>%
    filter(cancer == "COAD") %>%
    filter_depmap_by_allele_count(min = 3)
pryr::object_size(depmap_modelling_df)

################################################################################
################################################################################

set.seed(0)

GRAPHS_DIR <- "10_10_linear-modeling-syn-let"
# reset_graphs_directory(GRAPHS_DIR)




##### Add this to "lib/..."
# Center and standardize a numeric value
scale_numeric <- function(x, na.rm = FALSE) {
    (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
}



#### ---- Linear model of gene effect regressed on RNA expr. ---- ####

# Fit a linear model of gene effect on RNA expression.
rna_linear_model <- function(data) {
    if (all(is.na(data$rna_expression_std))) {
        return(NA)
    }

    stop("How to incorporate if the gene is mutated or not?")

    lm(gene_effect ~ 1 + rna_expression_std, data = data)
}


# Get the p-value of the linear model.
get_rna_linear_model_pvalue <- function(fit) {
    if (all(is.na(fit))) {
        return(1)
    } else {
        return(glance(fit)$p.value[[1]])
    }
}


# Make some diagnostic plots for the linear model of gene effect and RNA.
rna_linear_model_diagnostics <- function(hugo_symbol, fit) {
    p1 <- ggplot(fit, aes(.fitted, .resid)) +
        geom_point() +
        stat_smooth(method = "loess") +
        geom_hline(yintercept = 0, color = "red", lty = 2) +
        labs(x = "Fitted values",
             y = "Residuals",
             title = "Residual vs Fitted Plot")

    p2 <- ggplot(fit, aes(qqnorm(.stdresid)[[1]], .stdresid)) +
        geom_point(na.rm = TRUE) +
        geom_abline(slope = 1, intercept = 0) +
        labs(x = "Theoretical Quantiles",
             y = "Standardized Residuals",
             title = "Normal Q-Q")

    p3 <- ggplot(fit, aes(.fitted, sqrt(abs(.stdresid)))) +
        geom_point(na.rm = TRUE) +
        stat_smooth(method = "loess", na.rm = TRUE) +
        labs(x = "Fitted Value",
             y = expression(sqrt("|Standardized residuals|")),
             title = "Scale-Location")

    p4 <- ggplot(fit, aes(seq_along(.cooksd), .cooksd)) +
        geom_bar(stat = "identity", position = "identity") +
        labs(x = "Obs. Number",
             y = "Cook's distance",
             title = "Cook's distance")

    p5 <- ggplot(fit, aes(.hat, .stdresid)) +
        geom_point(aes(size = .cooksd), na.rm = TRUE) +
        stat_smooth(method = "loess", na.rm = TRUE) +
        labs(x = "Leverage",
             y = "Standardized Residuals",
             title = "Residual vs Leverage Plot") +
        scale_size_continuous("Cook's Distance", range = c(1, 5))

    p6 <- ggplot(fit, aes(.hat, .cooksd)) +
        geom_point(na.rm = TRUE) +
        stat_smooth(method = "loess", na.rm = TRUE) +
        labs(x = "Leverage hii",
             y = "Cook's Distance",
             title = "Cook's dist vs Leverage hii/(1-hii)") +
        geom_abline(slope = seq(0, 3, 0.5), intercept = 0,
                    color = "gray", lty = 2)

    full_p <- wrap_plots(list(p1, p2, p3, p4, p5, p6), nrow = 2) +
        plot_annotation(title = hugo_symbol) &
        theme_bw(base_size = 8, base_family = "Arial")  &
        theme(legend.position = "bottom")


    # TODO: need to change this in real script.
    ggsave(
        filename = file.path("graphs", GRAPHS_DIR,
                             glue("{hugo_symbol}_rna-model-diagnostics.jpeg")),
        plot = full_p,
        width = 12,
        height = 8,
        units = "in"
    )

    return(fit)
}



#### ---- KRAS allele ANOVA ---- ####

# Fit an ANOVA on gene effect with KRAS alleles as predictors.
kras_allele_anova <- function(data) {
    aov(gene_effect ~ kras_allele, data = data)
}


# Get the p-value of the ANOVA.
get_kras_allele_anova_pvalue <- function(fit) {
    glance(fit)$p.value[[1]]
}


# Get model diagnostics for the ANOVA.
get_kras_allele_anova_stats <- function(fit) {
    as_tibble(sjstats::anova_stats(fit))
}


# Plot some diagnostic plots for the ANOVA.
kras_allele_anova_diagnostics <- function(hugo_symbol, fit) {
    p1 <- ggplot(fit, aes(.fitted, .resid)) +
        geom_point() +
        stat_smooth(method = "lm") +
        geom_hline(yintercept = 0, color = "red", lty = 2) +
        labs(x = "Fitted values",
             y = "Residuals",
             title = "Residual vs Fitted Plot")
    p2 <- ggplot(fit$model, aes(x = kras_allele, y = gene_effect)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) +
        geom_hline(yintercept = 0, color = "red", lty = 2)

    full_p <- (p2 | p1) +
        plot_annotation(title = hugo_symbol) &
        theme_bw(base_size = 8, base_family = "Arial")  &
        theme(legend.position = "bottom")


    # TODO: need to change this in real script.
    ggsave(
        filename = file.path("graphs", GRAPHS_DIR,
                             glue("{hugo_symbol}_allele-anova-diagnositics.jpeg")),
        plot = full_p,
        width = 8,
        height = 4,
        units = "in"
    )

    return(fit)
}



temp_model_int <- depmap_modelling_df %>%
    group_by(cancer, hugo_symbol) %>%
    mutate(rna_expression_std = scale_numeric(rna_expression, na.rm = TRUE)) %>%
    group_by(cancer, hugo_symbol) %>%
    nest() %>%
    ungroup() %>%
    slice(1:5) %>%
    mutate(
        rna_lm = map(data, rna_linear_model),
        rna_lm = map2(hugo_symbol, rna_lm, rna_linear_model_diagnostics),
        rna_pvalue = map_dbl(rna_lm, get_rna_linear_model_pvalue)
    ) %>%
    mutate(
        allele_aov = map(data, kras_allele_anova),
        allele_aov_stats = map(allele_aov, get_kras_allele_anova_stats),
        allele_aov_pvalue = map_dbl(allele_aov, get_kras_allele_anova_pvalue),
        allele_aov = map2(hugo_symbol, allele_aov, kras_allele_anova_diagnostics)
    )
