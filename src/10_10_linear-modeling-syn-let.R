# Modeling the effect of knocking out each gene by the KRAS allele of the cell.
# Account for if the RNA expression or mutation status explains the effect.

set.seed(0)

GRAPHS_DIR <- "10_10_linear-modeling-syn-let"
reset_graph_directory(GRAPHS_DIR)



#### ---- Linear model of gene effect regressed on RNA expr. ---- ####

# Fit a linear model of gene effect on RNA expression.
rna_linear_model <- function(data) {
    if (all(is.na(data$rna_expression_std))) {
        return(NA)
    }
    lm(gene_effect ~ 1 + rna_expression_std, data = data)
}


# Get the p-value of the linear model.
get_rna_linear_model_stats <- function(fit) {
    if (all(is.na(fit))) {
        return(tibble())
    }
    glance(fit) %>%
        janitor::clean_names()
}


# Make some diagnostic plots for the linear model of gene effect and RNA.
rna_linear_model_diagnosticplots <- function(hugo_symbol, fit) {
    if (all(is.na(fit))) { return(fit) }

    p1 <- ggplot(fit, aes(.fitted, .resid)) +
        geom_point() +
        stat_smooth(method = "loess") +
        geom_hline(yintercept = 0, color = "red", lty = 2) +
        labs(x = "Fitted values",
             y = "Residuals",
             title = "Residual vs Fitted Plot")

    p2 <- ggplot(fit, aes(qqnorm(.stdresid, plot.it = FALSE)[[1]], .stdresid)) +
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


    filename <- file.path("graphs", GRAPHS_DIR,
                          glue("{hugo_symbol}_rna-model-diagnostics.svg"))
    ggsave_wrapper(full_p, filename, width = 12, height = 8)

    return(fit)
}


#### ---- Test if mutation matters ---- ####

# Test if whether the gene is mutated or not is associated with gene effect.
# The test is only run if there are at least `min_alterations` cell lines with
# a mutation.
gene_is_mutated_test <- function(data, min_alterations = 3) {
    q1 <- n_distinct(data$is_mutated) == 2
    q2 <- sum(data$is_mutated) >= min_alterations
    if (q1 & q2) {
        return(wilcox.test(gene_effect ~ is_mutated, data = data))
    }
    return(NA)
}


# Get the p-value from the test on whether the gene being mutated is
# associated with the gene effect.
get_gene_is_mutated_pvalue <- function(fit) {
    if (all(is.na(fit))) {
        return(1)
    } else {
        return(glance(fit)$p.value[[1]])
    }
}


# Plot some diagnositics for if the gene being mutated is associated with
# the gene effect.
gene_is_mutated_diagnostics <- function(hugo_symbol, data) {
    p <- data %>%
        mutate(is_mutated = as.character(is_mutated),
               is_mutated = factor(is_mutated, levels = c("TRUE", "FALSE"))) %>%
        ggplot(aes(x = is_mutated, y = gene_effect)) +
        geom_boxplot() +
        geom_jitter(width = 0.2) +
        theme_bw(base_size = 8, base_family = "Arial") +
        labs(title = hugo_symbol)

    filename <- file.path("graphs", GRAPHS_DIR,
                          glue("{hugo_symbol}_is-mutated-diagnostics.svg"))
    ggsave_wrapper(p, filename, "small")

    return(data)
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


    filename <- file.path("graphs", GRAPHS_DIR,
                          glue("{hugo_symbol}_allele-anova-diagnositics.svg"))
    ggsave_wrapper(full_p, filename, "wide")

    return(fit)
}



#### ---- One-vs.-All pairwise comparisons ---- ####

# Test the gene effect between each allele and the rest.
one_vs_all_pariwise_test <- function(data, padjust_method = "BH") {

    run_t_test <- function(allele) {
        d <- data %>% mutate(is_allele = kras_allele == !!allele)
        t.test(gene_effect ~ is_allele, data = d) %>%
            glance() %>%
            rename(estimate_other = estimate1,
                   estimate_allele = estimate2)
    }

    tibble(allele = sort(unique(data$kras_allele)),
           test = map(allele, run_t_test)) %>%
        unnest(test) %>%
        janitor::clean_names() %>%
        mutate(adj_p_value = p.adjust(p_value, method = padjust_method))
}



#### ---- Pairwise allele comparisons ---- ####

# Pairwise comparisons of gene effect between the KRAS alleles.
one_vs_one_pairwise_test <- function(data, padjust_method = "BH") {
    get_mean_value <- function(a) {
        data %>%
            filter(kras_allele == !!a) %>%
            pull(gene_effect) %>%
            unlist() %>%
            mean()
    }

    pairwise.wilcox.test(x = data$gene_effect,
                         g = data$kras_allele,
                         p.adjust.method = "none") %>%
        tidy() %>%
        janitor::clean_names() %>%
        mutate(adj_p_value = p.adjust(p_value, method = padjust_method),
               estimate1 = map_dbl(group1, get_mean_value),
               estimate2 = map_dbl(group2, get_mean_value),
               estimate = estimate1 - estimate2)
}



#### ---- Run the data throgh the workflow ---- ####

# Remove rows where the gene is deleted in the cell line.
# The entire gene is removed if this reduces the number of cell lines per cancer
# below 3 cell lines.

cache("depmap_model_workflow_res",
      depends = "depmap_modelling_df",
{
    depmap_model_workflow_res <- depmap_modelling_df %>%
        filter_depmap_by_allele_count() %>%
        group_by(cancer) %>%
        filter(n_distinct(kras_allele) >= 3) %>%
        filter(!is_deleted) %>%
        add_count(cancer, hugo_symbol, kras_allele) %>%
        group_by(cancer, hugo_symbol) %>%
        filter(all(n >= 3)) %>%
        select(-n) %>%
        ungroup()

    # Scale RNA expression and nest data per gene per cancer.
    depmap_model_workflow_res <- depmap_model_workflow_res %>%
        group_by(cancer, hugo_symbol) %>%
        mutate(rna_expression_std = scale_numeric(rna_expression,
                                                  na.rm = TRUE)) %>%
        group_by(cancer, hugo_symbol) %>%
        nest() %>%
        ungroup()

    # Run data through the modeling workflow,
    depmap_model_workflow_res <- depmap_model_workflow_res %>%
        mutate(cancer_hugo = paste(cancer, "-", hugo_symbol)) %>%
        mutate(
            rna_lm = map(data, rna_linear_model),
            rna_lm_stats = map(rna_lm, get_rna_linear_model_stats)
        ) %>%
        mutate(
            ismut_fit = map(data, gene_is_mutated_test),
            ismut_pvalue = map_dbl(ismut_fit, get_gene_is_mutated_pvalue)
        ) %>%
        mutate(
            allele_aov = map(data, kras_allele_anova),
            allele_aov_stats = map(allele_aov, get_kras_allele_anova_stats),
            allele_aov_pvalue = map_dbl(allele_aov,
                                        get_kras_allele_anova_pvalue)
        ) %>%
        mutate(
            ova_pairs = map(data, one_vs_all_pariwise_test),
            ova_any_sig = map_lgl(ova_pairs, ~ any(.x$adj_p_value < 0.05))
        ) %>%
        mutate(
            ovo_pairs = map(data, one_vs_one_pairwise_test),
            ovo_any_sig = map_lgl(ovo_pairs, ~ any(.x$adj_p_value < 0.05))
        )
})


#### ---- Inspection of genes that fit RNA linear model ---- ####

rna_lm_plot1 <- depmap_model_workflow_res %>%
    unnest(rna_lm_stats,  names_sep = "_") %>%
    mutate(is_sig = as.character(rna_lm_stats_p_value < 0.01)) %>%
    ggplot(aes(x = rna_lm_stats_r_squared, y = -log10(rna_lm_stats_p_value))) +
    facet_wrap(~ cancer, nrow = 1, scales = "free") +
    geom_point(aes(size = is_sig, color = is_sig)) +
    geom_hline(yintercept = -log10(0.05), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.01), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.001), lty = 2, size = 0.4) +
    scale_size_manual(values = c(0.05, 0.3)) +
    scale_color_manual(values = c("grey70", "grey25")) +
    theme_bw(base_size = 8, base_family = "Arial")
ggsave_wrapper(rna_lm_plot1,
               plot_path(GRAPHS_DIR, "rna_lm_plot1.svg"),
               "wide")


rna_lm_plot2 <- depmap_model_workflow_res %>%
    unnest(rna_lm_stats,  names_sep = "_") %>%
    mutate(
        is_sig = as.character(rna_lm_stats_p_value < 0.01),
        cancer_hugo = fct_reorder(cancer_hugo, -log10(rna_lm_stats_p_value))
    ) %>%
    ggplot(aes(x = cancer_hugo, y = -log10(rna_lm_stats_p_value))) +
    facet_wrap(~ cancer, nrow = 1, scales = "free") +
    geom_point(aes(size = is_sig, color = is_sig)) +
    geom_hline(yintercept = -log10(0.05), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.01), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.001), lty = 2, size = 0.4) +
    scale_size_manual(values = c(0.05, 0.3)) +
    scale_color_manual(values = c("grey70", "grey25")) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank()
    )
ggsave_wrapper(rna_lm_plot2,
               plot_path(GRAPHS_DIR, "rna_lm_plot2.svg"),
               "wide")


#### ---- Overview plots of genes with diff. dependency by KRAS ---- ####

# One vs. all volcano plot
ova_volcano_plot <- depmap_model_workflow_res %>%
    unnest(rna_lm_stats,  names_sep = "_") %>%
    filter(rna_lm_stats_p_value > 0.05) %>%
    filter(ismut_pvalue > 0.05) %>%
    select(cancer, hugo_symbol, cancer_hugo, ova_pairs) %>%
    unnest(ova_pairs) %>%
    filter(allele != "WT") %>%
    mutate(label = ifelse(-log10(adj_p_value) >= 3, hugo_symbol, NA)) %>%
    ggplot(aes(x = estimate, y = -log10(adj_p_value))) +
    facet_wrap(~ cancer, nrow = 1, scales = "free") +
    geom_point(aes(color = allele), size = 0.06, alpha = 0.5) +
    ggrepel::geom_text_repel(aes(label = label),
                             size = 0.8, force = 0.1, family = "Arial",
                             segment.size = 0.2) +
    scale_color_manual(values = short_allele_pal) +
    geom_hline(yintercept = -log10(0.05), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.01), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.001), lty = 2, size = 0.4) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        axis.ticks = element_blank()
    )
ggsave_wrapper(ova_volcano_plot,
               plot_path(GRAPHS_DIR, "ova_volcano_plot.svg"),
               "wide")


ovo_volcano_plot <- depmap_model_workflow_res %>%
    unnest(rna_lm_stats,  names_sep = "_") %>%
    filter(rna_lm_stats_p_value > 0.05) %>%
    filter(ismut_pvalue > 0.05) %>%
    select(cancer, hugo_symbol, cancer_hugo, ovo_pairs) %>%
    unnest(ovo_pairs) %>%
    filter(group1 != "WT" & group2 != "WT") %>%
    mutate(label = ifelse(-log10(adj_p_value) >= 3, hugo_symbol, NA),
           color = ifelse(adj_p_value < 0.01, "black", "grey70")) %>%
    ggplot(aes(x = estimate, y = -log10(adj_p_value))) +
    facet_grid(~ cancer, scales = "free") +
    geom_point(aes(color = color), size = 0.06) +
    ggrepel::geom_text_repel(aes(label = label),
                             size = 0.8, force = 0.1, family = "Arial",
                             segment.size = 0.2) +
    scale_color_identity() +
    geom_hline(yintercept = -log10(0.05), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.01), lty = 2, size = 0.4) +
    geom_hline(yintercept = -log10(0.001), lty = 2, size = 0.4) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        axis.ticks = element_blank()
    )
ggsave_wrapper(ovo_volcano_plot,
               plot_path(GRAPHS_DIR, "ovo_volcano_plot.svg"),
               "wide")
