
# Make a tibble to match names of genes.

cache("gene_name_tib", {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    gene_name_tib <- full_join(
            as_tibble(as.data.frame(org.Hs.egENSEMBL)),
            as_tibble(as.data.frame(org.Hs.egSYMBOL)),
            by = "gene_id"
        ) %>%
        dplyr::select(-gene_id) %>%
        filter(!is.na(ensembl_id) & !is.na(symbol))
    return(gene_name_tib)
})
