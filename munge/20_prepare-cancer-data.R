#################################
## Prepare GM data for my uses ##
#################################

# use ~18GB of RAM

mutation_type_regex <- "frame|missense|nonsense|splice|nonstop|start|targeted"

cache("cancer_muts_df",
{
    info(logger, "Preparing cancer data from GM.")

    get_data_path <- function(x) file.path("data", "cancer-data", x)

    # cancer data
    cancer_data <- get_data_path("wxgs_data_paad_luad_coadread_skcm_mm_Jun2019.rds") %>%
        readRDS()

    hypermutants <- get_data_path("hypermut_Aug2019.rds") %>%
        readRDS() %>%
        unlist() %>%
        unique()

    mutation_type_regex <- "frame|missense|nonsense|splice|nonstop|start|targeted"

    cancer_data_frames <- c()
    for (CANCER in names(cancer_data)) {
        cancer_name <- str_to_upper(ifelse(CANCER == "coadread", "coad", CANCER))

        df <- cancer_data[[CANCER]]$mutations$data %>%
            as_tibble() %>%
            mutate(cancer = !! cancer_name)
        cancer_data_frames <- c(cancer_data_frames, list(df))
    }
    rm(cancer_data)

    cancer_muts_df <- bind_rows(cancer_data_frames) %>%
        select(-entrez_gene_id, -tumor_type) %>%
        dplyr::rename(gene = "gene_symbol", dataset = "genetic_profile_id") %>%
        mutate(mutation_type = str_to_lower(mutation_type),
               is_hypermutant = case_id %in% hypermutants)


#### ---- Assigning KRAS mutation ---- ####

    ras_hotspots <- c("12", "13", "61", "146", "170", "117")
    ras_mutant_tib <- cancer_muts_df %>%
        filter(
            gene == "KRAS" &
            amino_position %in% !!ras_hotspots &
            str_detect(mutation_type, !!mutation_type_regex)
        ) %>%
        dplyr::rename(ras = "gene") %>%
        mutate(ras_allele = paste0(ras, "_", amino_acid_change)) %>%
        group_by(cancer, dataset, case_id) %>%
        filter(VAF == max(VAF) | is.na(VAF)) %>%
        ungroup() %>%
        unique()

    double_kras_mutants <- ras_mutant_tib %>%
        group_by(case_id) %>%
        filter(n() > 1) %>%
        ungroup() %>%
        jhcutils::u_pull(case_id) %T>%
        saveRDS(get_data_path("double_kras_mutants.rds"))

    ras_mutant_tib %<>%
        filter(!(case_id %in% !!double_kras_mutants))

    cancer_muts_df %<>%
        filter(!(case_id %in% !!double_kras_mutants))

    # should be 0 => no double mutants!
    num_ras_double_mutants_left <- ras_mutant_tib %>%
        group_by(case_id) %>%
        filter(n() > 1) %>%
        nrow()
    stopifnot(num_ras_double_mutants_left == 0)
    # should be 0

    # cache
    info(logger, "Caching RAS mutant table for cancer data.")
    cache("ras_mutant_tib", { return(ras_mutant_tib) })

    # trim to only required cells for join to `cancer_muts_df`
    ras_mutant_tib %<>%
        select(cancer, dataset, case_id, ras, ras_allele)

    cancer_muts_df %<>%
        left_join(ras_mutant_tib, by = c("cancer", "dataset", "case_id")) %>%
        mutate(ras = ifelse(is.na(ras), "WT", ras),
               ras_allele = ifelse(is.na(ras_allele), "WT", ras_allele))

    return(cancer_muts_df)
})


cache("cancer_coding_muts_df",
      depends = "cancer_coding_muts",
{
    # save only non-silent mutations
    cancer_coding_muts_df <- cancer_muts_df %>%
        filter(str_detect(mutation_type, mutation_type_regex))
    return(cancer_coding_muts_df)
})
