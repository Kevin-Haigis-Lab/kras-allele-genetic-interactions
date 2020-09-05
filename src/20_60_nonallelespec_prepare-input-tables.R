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
write_mutation_tables <- function(cancer, data, save_dir, ...) {
    file_name <- make_file_name(cancer, save_dir)
    data %>%
        select(sampleid, gene) %>%
        unique() %>%
        write_tsv(path = file_name)
    invisible(cancer)
}

# create the file name to write to
make_file_name <- function(cancer, dir) {
    file.path(
        dir,
        paste(cancer, "mutations.tsv", sep = "_")
    )
}


input_directory <- file.path("data", "rc-test-nonallelespec", "input")
dir.create(input_directory, showWarnings = FALSE)


# write out the mutation tables for each cancer
a <- cancer_coding_muts_df %>%
    filter(!is_hypermutant) %>%
    filter(cancer != "SKCM") %>%
    select(cancer, sampleid = tumor_sample_barcode, gene = hugo_symbol) %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(write_mutation_tables, save_dir = input_directory)

rm(a)
