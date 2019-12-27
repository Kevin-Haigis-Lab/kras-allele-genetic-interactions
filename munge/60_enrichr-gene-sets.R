# Prepare the gene sets used by Enrichr downloaded from their library:
#   https://amp.pharm.mssm.edu/Enrichr/#stats

# Extract the term from the line of the gene set file.
get_terms <- function(l) {
    unlist(str_split(l, "\t"))[[1]]
}

# Get the genes for a gene set from the line of the gene set file.
get_genes <- function(l) {
    a <- unlist(str_split(l, "\t"))[-1]
    a[a != ""]
}

# Parse a gene set file from Enrichr's library.
parse_geneset_file <- function(path) {
    lines <- readLines(path)
    terms <- purrr::map_chr(lines, get_terms)
    genes <- purrr::map(lines, get_genes)
    invisible(tibble(term = terms, genes = genes))
}


ProjectTemplate::cache("enrichr_genesets",
{
    enrichr_db_dir <- file.path("data", "enrichr-gene-sets")
    enrichr_geneset_paths <- list.files(enrichr_db_dir,
                                        full.names = TRUE,
                                        pattern = "txt$")


    enrichr_genesets <- tibble(
            path = enrichr_geneset_paths,
            datasource = file_sans_ext(basename(enrichr_geneset_paths))
        ) %>%
        mutate(data = map(path, parse_geneset_file)) %>%
        unnest(data) %>%
        mutate(std_term = map2_chr(term, datasource, mod_term_for_datasource),
               std_term = standardize_enricher_terms(std_term)) %>%
        select(-path)

    return(enrichr_genesets)
})

