
# Look for functional enrichment in the genes with detectable genetic
# interactions with KRAS alleles

library(enrichR)

dbs <- as_tibble(listEnrichrDbs())
# dbs$libraryName %>% sort()

enrichr_dbs <- c(
    "BioCarta_2016",
    "GO_Biological_Process_2018",
    "KEA_2015",
    "KEGG_2019_Human",
    "LINCS_L1000_Kinase_Perturbations_down",
    "LINCS_L1000_Kinase_Perturbations_up",
    "Panther_2016",
    "PPI_Hub_Proteins",
    "Reactome_2016",
    "Transcription_Factor_PPIs",
    "WikiPathways_2019_Human"
)


enrichr_wrapper <- function(gene_list) {
    if (n_distinct(gene_list) == 0) { return(NULL) }

    res <- enrichr(gene_list, enrichr_dbs) %>%
        enframe(name = "datasource", value = "data") %>%
        mutate(size_of_res = purrr::map_int(data, ~ nrow(.x))) %>%
        filter(size_of_res > 0) %>%
        select(-size_of_res) %>%
        unnest() %>%
        janitor::clean_names()

    return(res)
}
enrichr_wrapper <- memoise::memoise(enrichr_wrapper)

cache("enrichr_tib", depends = "genetic_interaction_df",
{
    enrichr_tib <- genetic_interaction_df %>%
        group_by(cancer, allele) %>%
        summarise(gene_list = list(hugo_symbol)) %>%
        ungroup() %>%
        mutate(enrichr_res = purrr::map(gene_list, enrichr_wrapper))
    return(enrichr_tib)
})



# Parse the "overlap" column to get number of genes in the gene set
get_enrichr_overlap_int <- function(overlap) {
    as.integer(str_split_fixed(overlap, "/", 2)[, 1])
}


# Write out the results from Enrichr.
#  pval: cut-off for adjusted p-value
#  min_overlap: minimum number of genes in the gene set
write_enrichr_results <- function(cancer, allele, gene_list, enrichr_res,
                                  pval = 0.05, min_overlap = -1) {
    xlsx_save_path <- file.path(
        "tables",
        "20_45_fxnal-enrich-genetic-interactions",
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




# a regular expression to use to remove uninteresting gene sets
uniteresting_terms_regex <- c(
      "cardiac", "heart", "muscle", "cardio",
      "disease", "chronic", "cancer",
      "basal", "axon", "amoeb", "glioma",
      "Hepatitis", "infection", "virus", "viral",
      "circulatory", "nervous", "neuro",
      "alcohol", "carcinoma", "microglia", "melanoma", "thyroid",
      "blastoma", "autism", "oocyte", "depression"
) %>%
    paste0(collapse = "|") %>%
    regex(ignore_case = TRUE, dotall = TRUE)


# Make a dot-plot of the functions from a data source enriched for a cancer.
dotplot_top_functions <- function(cancer, datasource, data) {
    mod_data <- data %>%
        filter(!str_detect(term, !!uniteresting_terms_regex)) %>%
        mutate(
            term = str_remove_all(term, "_Homo.sapiens+.*$") %>%
                str_remove_all("_96h.*$") %>%
                str_remove_all("[:space:]WP.*$") %>%
                str_remove_all("[:space:]\\(GO.*$") %>%
                str_replace_all("_", " ") %>%
                str_wrap(40)
        ) %>%
        complete(allele, term) %>%
        mutate(
            adjusted_p_value = ifelse(is.na(adjusted_p_value), 1, adjusted_p_value),
            n_overlap = ifelse(is.na(n_overlap), 0, n_overlap)
        )

    if (nrow(mod_data) == 0) { return() }

    min_size <- ifelse(min(mod_data$n_overlap) == 0,  0,  1 )
    min_alpha <- ifelse(min(mod_data$n_overlap) == 0,  0,  0.1)

    p <- mod_data %>%
        ggplot(
            aes(x = allele, y = term)
        ) +
        geom_point(
            aes(size = -log10(adjusted_p_value),
                alpha = n_overlap),
            color = "blue"
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
    save_path <- plot_path("20_45_fxnal-enrich-genetic-interactions",
                           glue("enrichr_{cancer}_{datasource}.svg"))
    ggsave_wrapper(p, save_path, "large")
}


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
    group_by(cancer, datasource) %>%
    nest() %>%
    purrr::pwalk(dotplot_top_functions)
