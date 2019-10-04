

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


lm_on_rna <- function(data, ...) {
    data_mod <- lm1_prepare_data(data)
    fit <- lm(gene_effect ~ rna_scaled, data = data_mod)
    return(fit)
}


anova_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue < 0.01) { return(NA) }

    data_mod <- lm1_prepare_data(data)
    res <- aov(gene_effect ~ allele, data = data_mod)
    return(res)
}


kruskal_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue < 0.01) { return(NA) }

    data_mod <- lm1_prepare_data(data)
    res <- kruskal.test(gene_effect ~ allele, data = data_mod)
    return(res)
}


pairwise_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue < 0.01) { return(NA) }

    data_mod <- lm1_prepare_data(data)
    res <- pairwise.t.test(data_mod$gene_effect,
                           data_mod$allele,
                           p.adjust.method = "BH")
    return(res)
}

model1_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    ungroup() %>%
    slice(1:1000) %>% # TEST
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )


plot_pairwise_test_results <- function(hugo_symbol, cancer, data, allele_aov, allele_pairwise, ...) {
    if (all(is.na(allele_aov)) | all(is.na(allele_pairwise))) { return() }

    if (!any(tidy(allele_pairwise) < 0.10)) { return() }

    stat_tib <- compare_means(
        gene_effect ~ allele, data = data,
        method = "t.test", p.adjust.method = "BH"
    ) %>%
        filter(p.adj < 0.10)

    stat_bar_height <- 0.08
    stat_bar_y_positions <- c(max(data$gene_effect) + stat_bar_height)
    for (i in seq(1, nrow(stat_tib))) {
        stat_bar_y_positions <- c(
            stat_bar_y_positions,
            stat_bar_y_positions[(i - 1)] + stat_bar_height
        )
    }
    stat_tib$y.position <- stat_bar_y_positions

    ggboxplot(
            data,
            x = "allele",
            y = "gene_effect",
            color = "allele",
            add = "jitter"
        ) +
        stat_pvalue_manual(stat_tib, label = "p.adj") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = glue("Cancer: {cancer}, target gene: {hugo_symbol}"),
            y = "depletion effect"
        )

    browser()
}


model1_tib %>%
    pwalk(plot_pairwise_test_results)
