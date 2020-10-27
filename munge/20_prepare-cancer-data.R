#################################
## Prepare GM data for my uses ##
#################################

mutation_type_regex <- "frame|missense|nonsense|splice|nonstop|start|targeted"

get_data_path <- function(x) file.path("data", "cancer-data", x)


get_mutation_dataframes_out_of_gm_structures <- function(file_path) {
  gm_data <- readRDS(file_path)

  dfs <- c()
  for (CANCER in names(gm_data)) {
    cancer_name <- str_to_upper(ifelse(CANCER == "coadread", "coad", CANCER))

    df <- gm_data[[CANCER]]$mutations$data %>%
      as_tibble() %>%
      mutate(cancer = !!cancer_name)
    dfs <- c(dfs, list(df))
  }

  return(dfs)
}


initial_preparation_of_gm_cancer_data <- function(dfs) {
  hypermutants <- get_data_path("hypermut_Aug2019.rds") %>%
    readRDS() %>%
    unlist() %>%
    unique()

  bind_rows(dfs) %>%
    select(-entrez_gene_id, -tumor_type) %>%
    dplyr::rename(
      hugo_symbol = "gene_symbol",
      dataset = "genetic_profile_id",
      tumor_sample_barcode = "case_id"
    ) %>%
    mutate(
      mutation_type = str_to_lower(mutation_type),
      mutation_type_hr = unlist(mapping_mutation_types_to_human_readable[mutation_type]),
      is_hypermutant = tumor_sample_barcode %in% hypermutants
    )
}


get_kras_mutant_tib <- function(df) {
  df %>%
    filter(
      hugo_symbol == "KRAS" &
        amino_position %in% !!kras_hotspot_codons$char &
        str_detect(mutation_type, !!mutation_type_regex)
    ) %>%
    dplyr::rename(ras = "hugo_symbol") %>%
    mutate(ras_allele = paste0(ras, "_", amino_acid_change)) %>%
    group_by(cancer, dataset, tumor_sample_barcode) %>%
    filter(VAF == max(VAF) | is.na(VAF)) %>%
    ungroup() %>%
    unique()
}


get_kras_double_mutants <- function(kras_df) {
  kras_df %>%
    group_by(tumor_sample_barcode) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    jhcutils::u_pull(tumor_sample_barcode)
}


check_for_kras_doublemutants <- function(kras_df) {
  # should be 0 --> no double mutants!
  num_ras_double_mutants_left <- ras_mutant_tib %>%
    group_by(tumor_sample_barcode) %>%
    filter(n() > 1) %>%
    nrow()
  if (num_ras_double_mutants_left != 0) {
    fatal(logger, "There are still double KRAS mutants.")
    return(FALSE)
  } else {
    info(logger, "There are no more KRAS double mutants.")
    return(TRUE)
  }
}


join_cancerdf_and_krasdf <- function(cancer_df, kras_df) {
  kras_df_mod <- kras_df %>%
    select(cancer, dataset, tumor_sample_barcode, ras, ras_allele)

  cancer_df_mod <- cancer_df %>%
    left_join(kras_df_mod, by = c("cancer", "dataset", "tumor_sample_barcode")) %>%
    mutate(
      ras = ifelse(is.na(ras), "WT", ras),
      ras_allele = ifelse(is.na(ras_allele), "WT", ras_allele)
    )
}


ProjectTemplate::cache("cancer_muts_df", {
  info(logger, "Preparing cancer data from GM.")

  # Read in and prepare the nested structure from GM.
  cancer_muts_df <- get_data_path(
    "wxgs_data_paad_luad_coadread_skcm_mm_Jun2019.rds"
  ) %>%
    get_mutation_dataframes_out_of_gm_structures() %>%
    initial_preparation_of_gm_cancer_data()

  # Assign KRAS mutations.
  ras_mutant_tib <- get_kras_mutant_tib(cancer_muts_df)
  double_kras_mutants <- get_kras_double_mutants(ras_mutant_tib)

  ras_mutant_tib %<>%
    filter(!(tumor_sample_barcode %in% !!double_kras_mutants))

  cancer_muts_df %<>%
    filter(!(tumor_sample_barcode %in% !!double_kras_mutants))

  check_for_kras_doublemutants(ras_mutant_tib)

  # trim to only required cells for join to `cancer_muts_df`
  cancer_muts_df <- join_cancerdf_and_krasdf(cancer_muts_df, ras_mutant_tib)

  return(cancer_muts_df)
})


ProjectTemplate::cache("cancer_coding_muts_df",
  depends = "cancer_muts_df",
  {
    # save only non-silent mutations
    cancer_coding_muts_df <- cancer_muts_df %>%
      filter(str_detect(mutation_type, mutation_type_regex))
    return(cancer_coding_muts_df)
  }
)




#### ---- Full cancer data set from GM ---- ####


ProjectTemplate::cache("cancer_full_muts_df", {
  info(logger, "Preparing cancer data from GM.")

  # Read in and prepare the nested structure from GM.
  cancer_full_muts_df <- get_data_path(
    "full_data_paad_luad_coadread_skcm_mm_Jun2019.rds"
  ) %>%
    get_mutation_dataframes_out_of_gm_structures() %>%
    initial_preparation_of_gm_cancer_data()

  # Assign KRAS mutations.
  ras_mutant_tib <- get_kras_mutant_tib(cancer_full_muts_df)
  double_kras_mutants <- get_kras_double_mutants(ras_mutant_tib)

  assign("double_kras_mutants", double_kras_mutants, envir = .GlobalEnv)
  ProjectTemplate::cache("double_kras_mutants")

  ras_mutant_tib %<>%
    filter(!(tumor_sample_barcode %in% !!double_kras_mutants))

  cancer_full_muts_df %<>%
    filter(!(tumor_sample_barcode %in% !!double_kras_mutants))

  check_for_kras_doublemutants(ras_mutant_tib)

  info(logger, "Caching RAS mutant table for cancer data.")
  assign("ras_mutant_tib", ras_mutant_tib, envir = .GlobalEnv)
  ProjectTemplate::cache("ras_mutant_tib")

  # trim to only required cells for join to `cancer_full_muts_df`
  cancer_full_muts_df <- join_cancerdf_and_krasdf(cancer_full_muts_df, ras_mutant_tib)

  return(cancer_full_muts_df)
})



ProjectTemplate::cache("cancer_full_coding_muts_df",
  depends = "cancer_full_muts_df",
  {
    # save only non-silent mutations
    cancer_full_coding_muts_df <- cancer_full_muts_df %>%
      filter(str_detect(mutation_type, mutation_type_regex))
    return(cancer_full_coding_muts_df)
  }
)
