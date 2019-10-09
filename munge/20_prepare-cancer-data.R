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
        dplyr::rename(
            hugo_symbol = "gene_symbol",
            dataset = "genetic_profile_id",
            tumor_sample_barcode = "case_id"
        ) %>%
        mutate(mutation_type = str_to_lower(mutation_type),
               is_hypermutant = tumor_sample_barcode %in% hypermutants)


#### ---- Assigning KRAS mutation ---- ####

    ras_mutant_tib <- cancer_muts_df %>%
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

    double_kras_mutants <- ras_mutant_tib %>%
        group_by(tumor_sample_barcode) %>%
        filter(n() > 1) %>%
        ungroup() %>%
        jhcutils::u_pull(tumor_sample_barcode)

    cache("double_kras_mutants", { return(double_kras_mutants) })

    ras_mutant_tib %<>%
        filter(!(tumor_sample_barcode %in% !!double_kras_mutants))

    cancer_muts_df %<>%
        filter(!(tumor_sample_barcode %in% !!double_kras_mutants))

    # should be 0 --> no double mutants!
    num_ras_double_mutants_left <- ras_mutant_tib %>%
        group_by(tumor_sample_barcode) %>%
        filter(n() > 1) %>%
        nrow()
    if (num_ras_double_mutants_left != 0) {
        fatal(logger, "There are still double KRAS mutants.")
    } else {
        info(logger, "There are no more KRAS double mutants.")
    }

    # cache
    info(logger, "Caching RAS mutant table for cancer data.")
    cache("ras_mutant_tib", { return(ras_mutant_tib) })

    # trim to only required cells for join to `cancer_muts_df`
    ras_mutant_tib %<>%
        select(cancer, dataset, tumor_sample_barcode, ras, ras_allele)

    cancer_muts_df %<>%
        left_join(ras_mutant_tib, by = c("cancer", "dataset", "tumor_sample_barcode")) %>%
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

