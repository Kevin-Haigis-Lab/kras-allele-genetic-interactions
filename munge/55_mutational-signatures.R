# Parse the mutational signatures results from GM.


ProjectTemplate::cache("mutational_signatures_df",
{
    # KRAS and tumor sample information.
    kras_muts <- cancer_muts_df %>%
        select(tumor_sample_barcode, dataset, target, cancer, is_hypermutant, ras_allele) %>%
        group_by(tumor_sample_barcode) %>%
        slice(1) %>%
        ungroup() %>%
        unique()

    # Build a tibble of signature contribution for each sample.
    mutational_signatures_df <- file.path(
        "data", "mutational-signatures",
        "signatures_cosmic_Doga_Aug2019.rds"
    ) %>%
        readRDS() %>%
        as.data.frame() %>%
        rownames_to_column("tumor_sample_barcode") %>%
        as_tibble() %>%
        pivot_longer(-tumor_sample_barcode,
                     names_to = "signature",
                     values_to = "contribution") %>%
        mutate(
            signature = str_remove(signature, "Signature\\."),
            signature = str_to_upper(signature)
        ) %>%
        left_join(kras_muts, by = "tumor_sample_barcode") %>%
        filter(!is_hypermutant & cancer != "SKCM") %>%
        left_join(signature_description_df, by = "signature")

    # Check no missing data for KRAS allele.
    if (any(is.na(mutational_signatures_df$ras_allele))) {
        stop("There are `NA` values for KRAS allele.")
    }

    # Check no missing signature descriptions.
    if (any(is.na(mutational_signatures_df$description))) {
        stop("Not all mut. sigs. have descriptions.")
    }

    log_rows(logger, mutational_signatures_df, "mutational_signatures_df")
    info(logger, "Caching `mutational_signatures_df`.")
    return(mutational_signatures_df)
})
