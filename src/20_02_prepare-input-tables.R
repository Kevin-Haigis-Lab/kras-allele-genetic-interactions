##################################################
## Create mutation tables for Snakemake-RC test ##
##################################################

# The Snakefile uses snakemake to run the RC test for mutual exclusivity and
# co-mutation on an input mutation table. This script makes the input mutation
# tables for each cancer.


# a function to write the mutation tables
# INPUT:
#    cancer: name of cancer
#    data: a tibble with columns "sampleid" and "gene"
#    min_mut_events: minimum number of mutations in the gene to be considered
# RETURNS:
#    writes a tsv file and (invisibly) returns the cancer name
write_mutation_tables <- function(cancer, data, ...) {
    file_name <- make_file_name(cancer)
    data %>%
        dplyr::select(sampleid, gene) %>%
        unique() %>%
        readr::write_tsv(path = file_name)
    invisible(cancer)
}

# create the file name to write to
make_file_name <- function(cancer) {
    file.path(
        "data", "rc-test", "input",
        paste(cancer, "mutations.tsv", sep = "_")
    )
}


# cancer data
cancer_data <- cancer_coding_muts_df %>%
    filter(!is_hypermutant)

# write out the mutation tables for each cancer
cancer_data %>%
    filter(cancer != "SKCM") %>%
    dplyr::mutate(gene = ifelse(gene == "KRAS", ras_allele, gene)) %>%
    dplyr::rename(sampleid = "case_id") %>%
    dplyr::select(cancer, sampleid, gene) %>%
    tidyr::nest(sampleid, gene) %>%
    purrr::pmap(write_mutation_tables)
