
# Look for functional enrichment in the genes with detectable genetic
# interactions with KRAS alleles

GRAPHS_DIR <- "20_45_fxnal-enrich-genetic-interactions"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


ProjectTemplate::cache("enrichr_tib",
                       depends = "genetic_interaction_df",
{
    enrichr_tib <- genetic_interaction_df %>%
        group_by(cancer, allele) %>%
        summarise(gene_list = list(hugo_symbol)) %>%
        ungroup() %>%
        mutate(enrichr_res = purrr::map(gene_list, enrichr_wrapper))
    return(enrichr_tib)
})


# Write out the results from Enrichr.
#  pval: cut-off for adjusted p-value
#  min_overlap: minimum number of genes in the gene set
write_enrichr_results <- function(cancer, allele, gene_list, enrichr_res,
                                  pval = 0.05, min_overlap = -1) {
    xlsx_save_path <- table_path(
        GRAPHS_DIR,
        glue("enrichr_results_{cancer}_{allele}.xlsx")
    )
    tsv_save_path <- str_replace(xlsx_save_path, "xlsx$", "tsv")
    res <- enrichr_res %>%
        mutate(n_genes = get_enrichr_overlap_int(overlap)) %>%
        filter(adjusted_p_value < !!pval & n_genes >= !!min_overlap) %>%
        select(-old_p_value, -old_adjusted_p_value) %>%
        arrange(adjusted_p_value, desc(n_genes))
    if (nrow(res) > 0) {
        res %T>%
            xlsx::write.xlsx(xlsx_save_path) %T>%
            write_tsv(tsv_save_path)
    }

}

purrr::pwalk(enrichr_tib, write_enrichr_results, pval = 0.2, min_overlap = 2)



# Make a dot-plot of the functions from a data source enriched for a cancer.
dotplot_top_functions <- function(cancer,
                                  datasource,
                                  data,
                                  reorder_terms = TRUE) {
    mod_data <- data %>%
        filter(!str_detect(term, !!uninteresting_enrichr_regex)) %>%
        mutate(
            term = str_wrap(term, 40),
            term = map2_chr(term, datasource, mod_term_for_datasource)
        ) %>%
        complete(allele, term) %>%
        mutate(
            adjusted_p_value = ifelse(is.na(adjusted_p_value), 1, adjusted_p_value),
            n_overlap = ifelse(is.na(n_overlap), 0, n_overlap)
        )

    if (nrow(mod_data) == 0) return()

    min_size <- ifelse(min(mod_data$n_overlap) == 0,  0,  1 )
    min_alpha <- ifelse(min(mod_data$n_overlap) == 0,  0,  0.1)

    if (reorder_terms) {
        term_levels <- mod_data %>%
            filter(!is.na(allele) & !is.na(adjusted_p_value) & n_overlap > 0) %>%
            group_by(term) %>%
            mutate(n_alleles = n_distinct(allele)) %>%
            ungroup() %>%
            arrange(desc(allele), n_alleles, -adjusted_p_value) %>%
            pull(term)

        mod_data$term <- factor(mod_data$term, levels = unique(term_levels))
    }

    # if (any(is.na(mod_data$allele))) {
    #     browser()
    # }

    # data %>%
    #     filter(is.na(allele)) %>%
    #     glimpse()

    p <- mod_data %>%
        ggplot(
            aes(x = allele, y = term)
        ) +
        geom_point(
            aes(size = -log10(adjusted_p_value),
                alpha = n_overlap),
            color = "dodgerblue"
        ) +
        scale_size_continuous(
            range = c(min_size, 8),
            guide = guide_legend(title.position = "top",
                                 title.hjust = 0.5,
                                 order = 10)
        ) +
        scale_alpha_continuous(
            range = c(min_alpha, 1),
            guide = guide_legend(title.position = "top",
                                 title.hjust = 0.5,
                                 order = 20)
        ) +
        theme_classic() +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            axis.title = element_blank()
        ) +
        labs(
            title = glue("{cancer} - data source: {datasource}"),
            alpha = "no. of genes",
            size = "-log10(adj. p-val)"
        )
    save_path <- plot_path(GRAPHS_DIR,
                           glue("enrichr_{cancer}_{datasource}.svg"))
    ggsave_wrapper(p, save_path, "large")
    invisible(p)
}


dotplot_selected_functions <- function(cancer, data) {
    p <- dotplot_top_functions(cancer, "SELECT", data)

    save_name <- glue("enrichr_{cancer}")
    if (cancer == "COAD") {
        saveRDS(p, get_fig_proto_path(save_name, 2))
    } else if (cancer == "LUAD") {
        saveRDS(p, get_fig_proto_path(save_name, 3))
    }
}

selected_enrichments <- tibble::tribble(
    ~cancer, ~datasource, ~term,
    "COAD", "WikiPathways_2019_Human", "Focal Adhesion",
    "COAD", "Panther_2016", "Hedgehog signaling pathway",
    "COAD", "Reactome_2016", "Beta-catenin phosphorylation cascade",
    "COAD", "Transcription_Factor_PPIs", "SMAD2",
    "COAD", "KEA_2015", "ATM",
    "COAD", 'PPI_Hub_Proteins', "YWHAZ",
    "COAD", "KEGG_2019_Human", "mTOR signaling pathway",
    "COAD", "KEGG_2019_Human", "Hippo signaling pathway",
    "COAD", "KEGG_2019_Human", "Cellular senescence",
    "COAD", "KEGG_2019_Human", "Apoptosis",
    "COAD", "KEGG_2019_Human", "PI3K-Akt signaling pathway",
    "COAD", "KEGG_2019_Human", "Wnt signaling pathway",
    "LUAD", "WikiPathways_2019_Human", "Focal Adhesion",
    "LUAD", "Panther_2016", "Wnt signaling pathway",
    "LUAD", "GO_Biological_Process_2018", "positive regulation of MAPK cascade",
    "LUAD", "Transcription_Factor_PPIs", "MYC",
    "LUAD", "WikiPathways_2019_Human", "RAC1/PAK1/p38/MMP2 Pathway",
    "LUAD", "WikiPathways_2019_Human", "Pathways Regulating Hippo Signaling",
    "LUAD", "PPI_Hub_Proteins", "PRKACA",
    "LUAD", "PPI_Hub_Proteins", "PRKCA",
    "LUAD", "PPI_Hub_Proteins", "CTNNB1",
    "LUAD", "KEGG_2019_Human", "PI3K-Akt signaling pathway",
    "LUAD", "KEGG_2019_Human", "Apelin signaling pathway",
    "LUAD", "BioCarta_2016", "Chromatin Remodeling by hSWI/SNF ATP-dependent Complexes",

)


enrichr_tib %>%
    select(-gene_list) %>%
    unnest(cols = enrichr_res) %>%
    mutate(n_overlap = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.05 & n_overlap >= 3) %>%
    mutate(overlap_genes = str_split(genes, ";")) %>%
    select(-genes) %>%
    group_by(term) %>%
    mutate(term_genes = list(unique(unlist(overlap_genes)))) %>%
    ungroup() %>%
    mutate(term = standardize_enricher_terms(term)) %>%
    group_by(cancer, datasource) %>%
    nest() %T>%
    purrr::pwalk(dotplot_top_functions) %>%
    unnest(data) %>%
    right_join(selected_enrichments,
               by = c("cancer", "datasource", "term")) %>%
    group_by(cancer) %>%
    nest() %>%
    purrr::pwalk(dotplot_selected_functions)
