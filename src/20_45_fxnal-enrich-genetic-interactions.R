
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

pwalk(enrichr_tib, write_enrichr_results, pval = 0.2, min_overlap = 2)


enrichr_significant_results <- enrichr_tib %>%
    select(-gene_list) %>%
    unnest(cols = enrichr_res) %>%
    mutate(n_overlap = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.05 & n_overlap >= 3) %>%
    mutate(overlap_genes = str_split(genes, ";")) %>%
    select(-genes) %>%
    group_by(term) %>%
    mutate(term_genes = list(unique(unlist(overlap_genes)))) %>%
    ungroup() %>%
    mutate(term = standardize_enricher_terms(term))


#### ---- Plot: Dot-plot per cancer and datasource ---- ####

enrichr_dotplot <- function(mod_data) {
    mod_data %>%
        ggplot(
            aes(x = allele, y = term)
        ) +
        geom_point(
            aes(alpha = -log10(adjusted_p_value),
                size = n_overlap),
            color = "dodgerblue"
        ) +
        theme_classic(base_size = 7, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            legend.title = element_markdown(),
            axis.title = element_blank(),
            strip.background = element_blank()
        ) +
        labs(
            size = "no. of genes",
            alpha = "-*log*<sub>10</sub>(adj. p-val)"
        )
}

get_term_order <- function(data) {
    data %>%
        filter(!is.na(allele) &
               !is.na(adjusted_p_value) &
               !is.na(n_overlap) &
               n_overlap > 0) %>%
        group_by(term) %>%
        mutate(n_alleles = n_distinct(allele)) %>%
        ungroup() %>%
        arrange(desc(allele), n_alleles, -adjusted_p_value) %>%
        pull(term)
}


# Make a dot-plot of the functions from a data source enriched for a cancer.
dotplot_top_functions <- function(cancer, datasource, data) {
    mod_data <- data %>%
        filter(!is.na(term))
    term_levels <- get_term_order(mod_data)
    mod_data$term <- factor(mod_data$term, levels = unique(term_levels))

    p <- enrichr_dotplot(mod_data) +
        scale_size_area(max_size = 5) +
        labs(title = glue("{cancer} - data source: {datasource}"))

    return(p)
}


save_plot <- function(plt, name, size = "large") {
    if (all(!is.na(plt))) {
        ggsave_wrapper(plt, plot_path(GRAPHS_DIR, name), size)
    }
    invisible(plt)
}


# Read in a data frame (tibble) with a manually-selected subset of
# enriched functions.
source(file.path("src", "20_44_select-enriched-functions.R"))

dotplot_single_datasource <- function(cancer, datasource, data) {
    return(tibble(
        cancer = cancer,
        datasource = datasource,
        data = list(data),
        plt = list(dotplot_top_functions(cancer,
                                         datasource,
                                         data[[1]]))
    ))
}


enrichr_significant_results %>%
    select(cancer, allele, datasource, term, n_overlap, adjusted_p_value) %>%
    filter(!str_detect(term, !!uninteresting_enrichr_regex)) %>%
    mutate(term = map2_chr(term, datasource, mod_term_for_datasource),
           term = str_wrap(term, 40)) %>%
    group_by(cancer, datasource) %>%
    nest() %>%
    ungroup() %>%
    pmap(dotplot_single_datasource) %>%
    bind_rows() %>%
    mutate(plt = map2(
        plt, paste0("enrichr_", cancer, "_", datasource, ".svg"), save_plot
    ))



#### ---- Plot: Dot-plot of select terms per cancer ---- ####


dotplot_selected_functions <- function(cancer, data) {
    p <- dotplot_top_functions(cancer, "SELECT", data)

    save_plot(p, as.character(glue("enrichr_{cancer}_SELECT.svg")))

    if (cancer %in% c("COAD", "LUAD", "PAAD")) {
        saveFigRds(p, as.character(glue("enrichr_{cancer}")))
    }
}


enrichr_significant_results %>%
    inner_join(selected_enrichments,
               by = c("cancer", "datasource", "term")) %>%
    mutate(term = map2_chr(term, datasource, mod_term_for_datasource),
           term = str_wrap(term, 40)) %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(dotplot_selected_functions)





#### ---- Plot: Dto-plot of select terms facet by cancer ---- ####

dotplot_selected_functions_faceted <- function(data) {
    mod_data <- data %>%
        mutate(term = paste0(term, "___", cancer))
    p <- dotplot_top_functions("FACET-CANCER", "SELECT", mod_data) +
        scale_y_discrete(labels = function(x) { str_remove(x, "___.*$") }) +
        facet_wrap(~ cancer, nrow = 1, scales = "free")
    save_plot(p, "enrichr_FACET_SELECT.svg")
    saveFigRds(p, "enrichr_all-cancers-faceted")
}


enrichr_significant_results %>%
    inner_join(selected_enrichments,
               by = c("cancer", "datasource", "term")) %>%
    mutate(term = map2_chr(term, datasource, mod_term_for_datasource),
           term = str_wrap(term, 30)) %>%
    dotplot_selected_functions_faceted()
