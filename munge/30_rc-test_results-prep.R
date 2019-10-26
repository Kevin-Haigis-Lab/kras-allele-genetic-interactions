
info(logger, "Combining all RC-test results into a single data frame.")

cache("rc_test_results",
      depends = "cancer_coding_muts_df",
{

    # get the other gene (not KRAS) in the `gene_set` column
    get_other_gene_from_genesets <- function(gss) {
        str_split_fixed(gss, " - ", 2) %>%
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
        file_parts <- basename(file_path) %>% str_split("_") %>% unlist()
        cancer <- file_parts[[1]]
        kras_allele <- file_parts[[2]] %>% str_replace("KRAS", "KRAS_")
        allele <- file_parts[[2]] %>% str_remove("KRAS")
        rc_test_type <- file_parts[[3]]
        df <- readRDS(file_path) %>%
            add_column(cancer = !!cancer,
                       kras_allele = !!kras_allele,
                       allele = !!allele,
                       rc_test_type = !!rc_test_type)
        return(df)
    }

    rc_test_results <- purrr::map(
            c("comutation", "exclusivity"),
            ~ list(list.files(file.path("data", "rc-test", "output", .x),
                         full.names = TRUE,
                         pattern = "rds"))
        ) %>%
        unlist() %>%
        purrr::map(read_rctest_output_file) %>%
        bind_rows() %>%
        mutate(hugo_symbol = get_other_gene_from_genesets(gene_sets)) %>%
        filter(!is.na(hugo_symbol)) %>%
        mutate(p_val = ifelse(p_val == 0, 1 / (2 * n_perms), p_val))


    num_samples_per_cancer_df <- cancer_coding_muts_df %>%
        filter(!is_hypermutant) %>%
        group_by(cancer) %>%
        summarise(
            num_samples_per_cancer = n_distinct(tumor_sample_barcode)
        )

    cancer_mut_counts <- cancer_coding_muts_df %>%
        filter(!is_hypermutant) %>%
        group_by(cancer, hugo_symbol) %>%
        summarise(
            num_mut_per_cancer = n_distinct(tumor_sample_barcode)
        ) %>%
        ungroup()


    num_samples_per_cancer_allele_df <- cancer_coding_muts_df %>%
        group_by(cancer, ras_allele) %>%
        summarise(
            num_samples_per_cancer_allele = n_distinct(tumor_sample_barcode)
        ) %>%
        ungroup()

    allele_mut_counts <- cancer_coding_muts_df %>%
        group_by(cancer, ras_allele, hugo_symbol) %>%
        summarise(
            num_mut_per_cancer_allele = n_distinct(tumor_sample_barcode)
        ) %>%
        ungroup()


    assign("cancer_mut_counts", cancer_mut_counts, envir = .GlobalEnv)
    cache("cancer_mut_counts")

    assign("allele_mut_counts", allele_mut_counts, envir = .GlobalEnv)
    cache("allele_mut_counts")


    rc_test_results <- rc_test_results %>%
        left_join(
            num_samples_per_cancer_df, by = "cancer"
        ) %>%
        left_join(
            cancer_mut_counts,
            by = c("cancer", "hugo_symbol")
        ) %>%
        mutate(num_mut_per_cancer = ifelse(
            is.na(num_mut_per_cancer), 0, num_mut_per_cancer)
        ) %>%
        left_join(
            num_samples_per_cancer_allele_df,
            by = c("cancer", "kras_allele" = "ras_allele")
        ) %>%
        left_join(allele_mut_counts,
                  by = c("cancer", "kras_allele" = "ras_allele", "hugo_symbol")
        ) %>%
        mutate(num_mut_per_cancer_allele = ifelse(
            is.na(num_mut_per_cancer_allele), 0, num_mut_per_cancer_allele
        ))


    return(rc_test_results)
})

