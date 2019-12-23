# Prepare the gene sets used by Enrichr downloaded from their library:
#   https://amp.pharm.mssm.edu/Enrichr/#stats


get_terms <- function(l) {
    unlist(str_split(l, "\t"))[[1]]
}

get_genes <- function(l) {
    a <- unlist(str_split(l, "\t"))[-1]
    a[a != ""]
}

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
        select(-path)

    return(enrichr_genesets)
})

