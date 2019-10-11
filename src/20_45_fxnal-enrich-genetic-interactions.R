
# Look for functional enrichment in the genes with detectable genetic
# interactions with KRAS alleles

library(enrichR)

dbs <- as_tibble(listEnrichrDbs())
dbs$libraryName %>% sort()

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


enrichr_tib <- genetic_interaction_df %>%
    group_by(cancer, allele) %>%
    summarise(gene_list = list(hugo_symbol)) %>%
    ungroup() %>%
    mutate(enrichr_res = purrr::map(gene_list, enrichr_wrapper))



# Parse the "overlap" column to get number of genes in the geneset
get_enrichr_overlap_int <- function(overlap) {
    as.integer(str_split_fixed(overlap, "/", 2)[, 1])
}


# Write out the results from Enrichr
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
