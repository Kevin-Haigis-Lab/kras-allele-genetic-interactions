
#### ---- Functional annotation of heatmap gene (row) clusters ---- ####

cluster_terms <- depmap_gene_clusters %>%
    group_by(cancer, gene_cls) %>%
    summarise(genes = list(hugo_symbol)) %>%
    ungroup() %>%
    mutate(enrichr_res = purrr::map(genes, enrichr_wrapper)) %>%
    select(-genes) %>%
    unnest(enrichr_res) %>%
    filter(!str_detect(term, !!uniteresting_enrichr_regex)) %>%
    mutate(n_genes = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.2 & n_genes > 2)

cache("cluster_terms", depends = "model1_tib")


# Write the enrichment results to file.
cluster_terms %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "gene-clusters-fxnal-enrichment.tsv"
    ))

cluster_terms %>%
    group_by(cancer, term) %>%
    filter(n_distinct(gene_cls) == 1) %>%
    group_by(cancer, datasource, gene_cls) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:10) %>%
    ungroup() %>%
    arrange(cancer, gene_cls, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "gene-clusters-fxnal-enrichment_top10.tsv"
    ))


common_term_regex <- c(
    "mitoch", "citric", "transcript", "translat"
) %>%
    paste0(collapse = "|") %>%
    regex(ignore_case = TRUE)

cluster_terms %>%
    filter(!str_detect(term, common_term_regex)) %>%
    group_by(cancer, term) %>%
    filter(n_distinct(gene_cls) == 1) %>%
    group_by(cancer, datasource, gene_cls) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:5) %>%
    ungroup() %>%
    arrange(cancer, gene_cls, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "gene-clusters-fxnal-enrichment_uncommon.tsv"
    ))


#### ---- Plot results ---- ####

# A side-ways bar-plot showing the adj. p-value and genes in enriched groups.
functional_enrichment_barplot <- function(cancer, data, genes_names_size = 3) {
    p <- data %>%
        ggplot(aes(x = term, y = -log10(adjusted_p_value))) +
        geom_col(aes(fill = -log10(adjusted_p_value)), alpha = 0.6) +
        geom_text(aes(label = genes),
                  color = "black",
                  y = 0, family = "Arial", hjust = 0, size = genes_names_size) +
        geom_text(aes(label = gene_cls),
                  hjust = 0.0, size = 3, family = "Arial", nudge_y = 0.01) +
        scale_fill_viridis_c() +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
        theme_bw(base_family = "Arial", base_size = 9) +
        theme(
            axis.title.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = glue("Functional enrichment in {cancer} synthetic lethals"),
            y = "-log10( adj. p-value )"
        ) +
        coord_flip()
    return(p)
}

# Make the side-ways bar-plot for the enriched functions of a gene-cluster
gene_cluster_functional_enrichment_barplot <- function(cancer,
                                                       data,
                                                       top_n_fxn = 10,
                                                       ...) {
    mod_data <- data %>%
        group_by(genes) %>%
        top_n(2, adjusted_p_value) %>%
        ungroup() %>%
        mutate(term = paste0(datasource, "-", term),
               term = str_wrap(term, 50),
               term = fct_reorder(term, -adjusted_p_value),
               genes = str_replace_all(genes, ";", ", ")) %>%
        arrange(adjusted_p_value, n_genes) %>%
        slice(seq(1, top_n_fxn))
    p <- functional_enrichment_barplot(cancer, mod_data)
    save_name <- plot_path("10_13_linear-modeling-syn-let_fxnal-enrichment",
                           glue("functional-enrichment_{cancer}.svg"))
    if (cancer == "LUAD") {
        ggsave_wrapper(p, save_name, "wide")
    } else {
        ggsave_wrapper(p, save_name, width = 6, height = 2)
    }
}

cluster_terms %>%
    filter(datasource != "LINCS_L1000_Kinase_Perturbations_down") %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(gene_cluster_functional_enrichment_barplot)



# Bar-plot of the gene effect scores for genes in enriched groups.
enriched_group_effect_barplot <- function(cancer, term, gene_cls, datasource,
                                          adjusted_p_value, odds_ratio, genes,
                                          ...) {
    p <- model_data %>%
        filter(cancer == !!cancer & hugo_symbol %in% enrichr_genes(genes)) %>%
        group_by(hugo_symbol, allele) %>%
        summarise(avg_gene_effect = mean(gene_effect),
                  sd_gene_effect = sd(gene_effect)) %>%
        ungroup() %>%
        mutate(bar_top_y = avg_gene_effect + sd_gene_effect,
               bar_bottom_y = avg_gene_effect - sd_gene_effect) %>%
        group_by(hugo_symbol) %>%
        mutate(index_value = mean(avg_gene_effect)) %>%
        ungroup() %>%
        mutate(hugo_symbol = fct_reorder(hugo_symbol, index_value)) %>%
        ggplot(aes(x = hugo_symbol, color = allele, fill = allele)) +
        geom_col(aes(y = avg_gene_effect), position = "dodge", alpha = 0.5) +
        geom_errorbar(aes(ymin = bar_bottom_y, ymax = bar_top_y),
                      position = position_dodge2(width = 0.25, padding = 0.75)) +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        scale_y_continuous(expand = expand_scale(mult = c(0.02, 0.02))) +
        theme_bw(base_family = "arial") +
        theme(
            axis.title.x = element_blank()
        ) +
        labs(
            title = glue("{cancer}: {term}"),
            y = "depletion effect"
        )

    save_term <- paste0(datasource, "_", term) %>%
        str_replace_sp() %>%
        str_replace_all("/", "-") %>%
        str_replace_all("\\(", "-") %>%
        str_remove_all("\\)|:")
    save_name <- plot_path(
        "10_13_linear-modeling-syn-let_fxnal-enrichment",
        glue("gene-effect-barplot_{cancer}_{save_term}.svg")
    )
    ggsave_wrapper(p, save_name, "wide")
}

cluster_terms %>%
    filter(datasource != "LINCS_L1000_Kinase_Perturbations_down") %>%
    pwalk(enriched_group_effect_barplot)



#### ---- Functional enrichment for each cancer (not gene clusters) ---- ####

cluster_terms_cancer <- depmap_gene_clusters %>%
    group_by(cancer) %>%
    summarise(genes = list(hugo_symbol)) %>%
    ungroup() %>%
    mutate(enrichr_res = purrr::map(genes, enrichr_wrapper)) %>%
    select(-genes) %>%
    unnest(enrichr_res) %>%
    filter(!str_detect(term, !!uniteresting_enrichr_regex)) %>%
    mutate(n_genes = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.2 & n_genes > 2)

cluster_terms_cancer %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "cancer-fxnal-enrichment.tsv"
    ))


cluster_terms_cancer %>%
    group_by(cancer, term) %>%
    group_by(cancer, datasource) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:10) %>%
    ungroup() %>%
    arrange(cancer, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "cancer-fxnal-enrichment_top10.tsv"
    ))

cluster_terms_cancer %>%
    filter(!str_detect(term, common_term_regex)) %>%
    group_by(cancer, datasource) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:10) %>%
    ungroup() %>%
    arrange(cancer, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "cancer-fxnal-enrichment_uncommon.tsv"
    ))


# Make the side-ways bar-plot for the enriched functions of a cancer
cancer_functional_enrichment_barplot <- function(cancer, data, ...) {
    mod_data <- data %>%
        group_by(genes) %>%
        top_n(2, adjusted_p_value) %>%
        ungroup() %>%
        mutate(term = paste(datasource, "-", term),
               term = str_wrap(term, 50),
               term = fct_reorder(term, -adjusted_p_value),
               genes = str_replace_all(genes, ";", ", ")) %>%
        arrange(adjusted_p_value, n_genes) %>%
        mutate(gene_cls = NA)
    p <- functional_enrichment_barplot(
        cancer, mod_data,
        genes_names_size = ifelse(cancer == "LUAD", 1.5, 3)
    )
    save_name <- plot_path("10_13_linear-modeling-syn-let_fxnal-enrichment",
                           glue("functional-enrichment_{cancer}_overall.svg"))
    if (cancer == "LUAD") {
        ggsave_wrapper(p, save_name, "large")
    } else {
        ggsave_wrapper(p, save_name, "wide")
    }

}

# make bar-plots
cluster_terms_cancer %>%
    filter(datasource != "LINCS_L1000_Kinase_Perturbations_down") %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(cancer_functional_enrichment_barplot)
