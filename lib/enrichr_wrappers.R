## Wrappers around the 'enrichR' pacakge for functional annotation

library(enrichR)

# A list of all databases available in Enrichr
all_enrichr_dbs <- as_tibble(listEnrichrDbs())


# Enrichr data bases used.
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

mapping_datasource_names <- list(
    "BioCarta_2016" = "BioCarta",
    "GO_Biological_Process_2018" = "GO_BP",
    "KEA_2015" = "KEA",
    "KEGG_2019_Human" = "KEGG",
    "LINCS_L1000_Kinase_Perturbations_down" = "LINCS_down",
    "LINCS_L1000_Kinase_Perturbations_up" = "LINCS_up",
    "Panther_2016" = "Panther",
    "PPI_Hub_Proteins" = "Hub",
    "Reactome_2016" = "Reactome",
    "Transcription_Factor_PPIs" = "TF",
    "WikiPathways_2019_Human" = "WikiPathways"
)

# Parse the "overlap" column to get number of genes in the gene set
get_enrichr_overlap_int <- function(overlap) {
    as.integer(str_split_fixed(overlap, "/", 2)[, 1])
}


# Get genes as a list.
enrichr_genes <- function(genes) {
    unlist(str_split(genes, ";"))
}


# A thin wrapper around Enrichr (memoised).
#   Supply a list of genes in `gene_list`.
#   Returns a tibble of the results.
#   The function is memoised to save time.
enrichr_wrapper <- function(gene_list) {
    if (n_distinct(gene_list) == 0) { return(NULL) }

    res <- enrichr(gene_list, enrichr_dbs) %>%
        enframe(name = "datasource", value = "data") %>%
        mutate(size_of_res = purrr::map_int(data, ~ nrow(.x))) %>%
        filter(size_of_res > 0) %>%
        select(-size_of_res) %>%
        unnest(cols = data) %>%
        janitor::clean_names()

    return(res)
}
enrichr_wrapper <- memoise::memoise(enrichr_wrapper)


# A regular expression to use to remove uninteresting gene sets
uninteresting_enrichr_regex <- c(
      "cardiac", "heart", "muscle", "cardio",
      "disease", "chronic", "cancer",
      "basal", "axon", "amoeb", "glioma",
      "Hepatitis", "infection", "virus", "viral",
      "circulatory", "nervous", "neuro",
      "alcohol", "carcinoma", "microglia", "melanoma", "thyroid",
      "blastoma", "autism", "oocyte", "depression", "spinal",
      "syndrome", "photodynamic", "caffeine", "HIV"
) %>%
    paste0(collapse = "|") %>%
    regex(ignore_case = TRUE, dotall = TRUE)
