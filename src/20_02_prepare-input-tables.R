##################################################
## Create mutation tables for Snakemake-RC test ##
##################################################

# The Snakefile uses snakemake to run the RC test for mutual exclusivity and
# co-mutation on an input mutation table. This script makes the input mutation
# tables for each cancer.


# Write the mutation tables.
write_mutation_tables <- function(cancer, data, save_dir, ...) {
  file_name <- make_file_name(cancer, save_dir)
  data %>%
    dplyr::select(sampleid, gene) %>%
    unique() %>%
    readr::write_tsv(path = file_name)
  invisible(cancer)
}

# Create the file name to write to.
make_file_name <- function(cancer, dir) {
  file.path(
    dir,
    paste(cancer, "mutations.tsv", sep = "_")
  )
}


input_directory <- here::here("data", "rc-test", "input")
if (!dir.exists(input_directory)) {
  message("Creating directory for input files.")
  dir.create(input_directory)
}


# write out the mutation tables for each cancer
a <- cancer_coding_muts_df %>%
  filter(!is_hypermutant) %>%
  filter(cancer != "SKCM") %>%
  dplyr::mutate(gene = ifelse(
    hugo_symbol == "KRAS", ras_allele, hugo_symbol
  )) %>%
  dplyr::rename(sampleid = "tumor_sample_barcode") %>%
  dplyr::select(cancer, sampleid, gene) %>%
  dplyr::group_by(cancer) %>%
  tidyr::nest() %>%
  purrr::pwalk(write_mutation_tables, save_dir = input_directory)

rm(a)
