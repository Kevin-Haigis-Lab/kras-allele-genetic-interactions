

essential_cutoffs <- model_data %>%
    filter(cancer != "MM") %>%
    select(cancer, cutoff_between_non_essential) %>%
    unique()


sort_ras_alleles <- function(alleles) {
    alleles <- unique(alleles)
    if (any(alleles == "WT")) {
        not_wt_alleles <- alleles[alleles != "WT"]
        return(c(
            sort(not_wt_alleles), "WT"
        ))
    } else {
        return(sort(alleles))
    }
}


pairwise_plot_wrapper <- function(cancer, hugo_symbol) {
    browser()

    df <- model_data %>%
        filter(cancer == !!cancer & hugo_symbol == !!hugo_symbol) %>%
        mutate(allele = factor(allele, levels = sort_ras_alleles(allele)))


    comparisons <- expand.grid(levels(df$allele), levels(df$allele))

    p <- ggboxplot(df
            x = "allele", y = "gene_effect",
            color = "allele", palette = short_allele_pal,
            add = "jitter"
        ) +
        stat_compare_means()
}


genes_to_plot <- model_data %>%
    right_join(depmap_gene_clusters, by = c("cancer", "hugo_symbol")) %>%
    group_by(cancer, hugo_symbol, allele) %>%
    summarise(avg_gene_effect = mean(gene_effect)) %>%
    ungroup() %>%
    left_join(essential_cutoffs, by = "cancer") %>%
    group_by(cancer, hugo_symbol) %>%
    filter(
        any(avg_gene_effect < cutoff_between_non_essential) &
        !(all(avg_gene_effect < cutoff_between_non_essential))
    ) %>%
    ungroup() %>%
    select(cancer, hugo_symbol) %>%
    unique() %>%
    pwalk(pairwise_plot_wrapper)