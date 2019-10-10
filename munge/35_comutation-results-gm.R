
info(logger, "Beginning munging of Fisher's comutation results.")

cache("fisher_comut_df", {
    fisher_comut_file_path <- file.path(
        "data", "fisher-comutation", "OR_fisher_multiTest_Aug2019.txt"
    )

    fisher_comut_df <- read_tsv(fisher_comut_file_path,
                                col_types = cols(),
                                progress = FALSE) %>%
        janitor::clean_names() %>%
        dplyr::rename(
            hugo_symbol = "gene",
            kras_allele = "allele",
            odds_ratio = "or",
            cancer = "tumor_type"
        ) %>%
        filter(
            cancer != "skcm" & str_detect(kras_allele, "KRAS")
        ) %>%
        mutate(
            kras_allele = jhcutils::str_replace_sp(kras_allele),
            cancer = str_to_upper(cancer),
            cancer = str_remove(cancer, "READ")
        ) %>%
        filter(kras_allele != "KRAS")

    return(fisher_comut_df)
})
