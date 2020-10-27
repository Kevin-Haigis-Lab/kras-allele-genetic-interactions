# Functions to help with processing RC test results.

# get the other gene (not KRAS) in the `gene_set` column
get_other_gene_from_genesets <- function(gss) {
  gss %>%
    str_split_fixed(" - ", 2) %>%
    as.data.frame() %>%
    set_names(c("V1", "V2")) %>%
    as_tibble() %>%
    mutate(keep_V1 = case_when(
      str_detect(V1, "KRAS") & str_detect(V2, "KRAS") ~ NA,
      str_detect(V1, "KRAS") ~ FALSE,
      str_detect(V2, "KRAS") ~ TRUE
    )) %>%
    mutate(other_gene = ifelse(keep_V1, V1, V2)) %>%
    pull(other_gene) %>%
    unlist()
}


# part the file name to get meta data
# read in the data frame and add the meta data as new columns
read_rctest_output_file <- function(file_path) {
  file_parts <- basename(file_path) %>%
    str_split("_") %>%
    unlist()
  cancer <- file_parts[[1]]
  kras_allele <- file_parts[[2]] %>% str_replace("KRAS", "KRAS_")
  allele <- file_parts[[2]] %>% str_remove("KRAS")
  rc_test_type <- file_parts[[3]]
  df <- readRDS(file_path) %>%
    add_column(
      cancer = !!cancer,
      kras_allele = !!kras_allele,
      allele = !!allele,
      rc_test_type = !!rc_test_type
    )
  return(df)
}


# First workflow for processing RC-test results.
initial_rc_results_processing <- function(file_list) {
  file_list %>%
    unlist() %>%
    purrr::map(read_rctest_output_file) %>%
    bind_rows() %>%
    mutate(hugo_symbol = get_other_gene_from_genesets(gene_sets)) %>%
    filter(!is.na(hugo_symbol)) %>%
    mutate(p_val = ifelse(p_val == 0, 1 / (2 * n_perms), p_val))
}


# Add some columns with relevant counts.
join_relevant_counts <- function(rc_res, muts_df) {
  num_samples_per_cancer_df <- muts_df %>%
    filter(!is_hypermutant) %>%
    group_by(cancer) %>%
    summarise(
      num_samples_per_cancer = n_distinct(tumor_sample_barcode)
    )

  cancer_mut_counts <- muts_df %>%
    filter(!is_hypermutant) %>%
    group_by(cancer, hugo_symbol) %>%
    summarise(
      num_mut_per_cancer = n_distinct(tumor_sample_barcode)
    ) %>%
    ungroup()


  num_samples_per_cancer_allele_df <- muts_df %>%
    group_by(cancer, ras_allele) %>%
    summarise(
      num_samples_per_cancer_allele = n_distinct(tumor_sample_barcode)
    ) %>%
    ungroup()

  allele_mut_counts <- muts_df %>%
    group_by(cancer, ras_allele, hugo_symbol) %>%
    summarise(
      num_mut_per_cancer_allele = n_distinct(tumor_sample_barcode)
    ) %>%
    ungroup()

  rc_res %>%
    left_join(
      num_samples_per_cancer_df,
      by = "cancer"
    ) %>%
    left_join(
      cancer_mut_counts,
      by = c("cancer", "hugo_symbol")
    ) %>%
    mutate(
      num_mut_per_cancer = ifelse(
        test = is.na(num_mut_per_cancer),
        yes = 0,
        no = num_mut_per_cancer
      )
    ) %>%
    left_join(
      num_samples_per_cancer_allele_df,
      by = c("cancer", "kras_allele" = "ras_allele")
    ) %>%
    left_join(allele_mut_counts,
      by = c("cancer", "kras_allele" = "ras_allele", "hugo_symbol")
    ) %>%
    mutate(
      num_mut_per_cancer_allele = ifelse(
        test = is.na(num_mut_per_cancer_allele),
        yes = 0,
        no = num_mut_per_cancer_allele
      )
    )
}
