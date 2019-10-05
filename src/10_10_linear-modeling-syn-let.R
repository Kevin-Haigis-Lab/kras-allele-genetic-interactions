
library(ggpubr)
library(glue)

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


rna_pvalue_is_significant <- function(pval, cutoff = 0.01) {
    if (is.na(pval)) { return(FALSE) }
    if (pval < cutoff) { return(TRUE) }
    return(FALSE)
}


anova_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    data_mod <- lm1_prepare_data(data)
    res <- aov(gene_effect ~ allele, data = data_mod)
    return(res)
}


kruskal_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    data_mod <- lm1_prepare_data(data)
    res <- kruskal.test(gene_effect ~ allele, data = data_mod)
    return(res)
}


pairwise_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    data_mod <- lm1_prepare_data(data)
    res <- pairwise.t.test(data_mod$gene_effect,
                           data_mod$allele,
                           p.adjust.method = "BH")
    return(res)
}

model1_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )


plot_pairwise_test_results <- function(hugo_symbol, cancer, data, allele_aov, allele_pairwise, ...) {
    if (all(is.na(allele_aov)) | all(is.na(allele_pairwise))) { return() }

    if (tidy(allele_aov)$p.value[[1]] >= 0.01) { return() }

    stat_tib <- compare_means(
        gene_effect ~ allele, data = data,
        method = "t.test", p.adjust.method = "BH"
    ) %>%
        filter(p.adj < 0.05)

    if (nrow(stat_tib) < 1) { return() }

    stat_bar_height <- 0.08
    stat_bar_y_positions <- c(max(data$gene_effect) + stat_bar_height)
    for (i in seq(1, nrow(stat_tib))) {
        stat_bar_y_positions <- c(
            stat_bar_y_positions,
            stat_bar_y_positions[(i - 1)] + stat_bar_height
        )
    }

    stat_tib$y.position <- stat_bar_y_positions

    p <- ggboxplot(
            data,
            x = "allele",
            y = "gene_effect",
            color = "allele",
            add = "jitter"
        ) +
        stat_pvalue_manual(stat_tib, label = "p.adj") +
        scale_color_manual(values = short_allele_pal) +
        theme(
            text = element_text(family = "arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = glue("Cancer: {cancer}, gene: {hugo_symbol}"),
            caption = "*bars indicate FDR-adjusted p-values < 0.05",
            y = "depletion effect"
        )

    plot_fname <- file.path(plot_save_dir, glue("{cancer}-{hugo_symbol}.svg"))
    ggsave_wrapper(p, plot_fname, size = "small")
}


# directory for save images
plot_save_dir <- plot_path("10_10_linear-modeling-syn-let_bosplots")
if (!dir.exists(plot_save_dir)) {
    info(logger, glue("Making directory for saving boxplots: {plot_save_dir}"))
    dir.create(plot_save_dir)
}

model1_tib %>%
    pwalk(plot_pairwise_test_results)
