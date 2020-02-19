# Prepare data of cancer prevalence.

DATA_DIR <- file.path("data", "cancer-rates")


# SOURCE: https://seer.cancer.gov/csr/1975_2016/browse_csr.php?sectionSEL=1&pageSEL=sect_01_table.04
ProjectTemplate::cache("seer_cancer_incidence_df", {
    usa_incidence_path <- file.path(
        DATA_DIR,
        "table1-4_incidence-death-survival-rates-all-races.csv"
    )

    seer_cancer_incidence_df <- read_csv(usa_incidence_path,
                                         skip = 3,
                                         n_max = 85) %>%
        janitor::clean_names() %>%
        select(site, contains("both_sexes"))

    colnames(seer_cancer_incidence_df) <- c(
        "site",
        "incidence_2012_2016",
        "us_mortality_2012_2016",
        "survival_percent_2009_2015"
    )

    cancer_names <- tibble::tribble(
                         ~site, ~cancer,
                   "all sites",   "all",
             "colon & rectum:",  "COAD",
             "lung & bronchus",  "LUAD",
                     "myeloma",    "MM",
                    "pancreas",  "PAAD",
        "melanoma of the skin",  "SKCM"
    )

    SCALE_FCTR <- 1e5

    seer_cancer_incidence_df %<>%
         mutate(
            site = str_to_lower(site),
            incidence_2012_2016 = as.numeric(incidence_2012_2016) * !!SCALE_FCTR,
            us_mortality_2012_2016 = as.numeric(us_mortality_2012_2016) * !!SCALE_FCTR,
            survival_percent_2009_2015 = as.numeric(survival_percent_2009_2015) * !!SCALE_FCTR
        ) %>%
        left_join(cancer_names, by = "site")

    return(seer_cancer_incidence_df)
})


# Preparation of ACS data.
# data_col_name: name to use for the data column (enquoted).
prepare_acs_data <- function(df, data_col_name) {
    data_col_name <- rlang::enquo(data_col_name)

    cancer_names <- tibble::tribble(
                       ~cancer_type, ~cancer,
        "all cancer types combined",   "all",
                       "colorectum",  "COAD",
                "lung and bronchus",  "LUAD",
                          "myeloma",    "MM",
                         "pancreas",  "PAAD",
             "melanoma of the skin",  "SKCM"
    )

    df %>%
        janitor::clean_names() %>%
        dplyr::rename(!!data_col_name := both_sexes_combined) %>%
        mutate(
            cancer_type = str_to_lower(cancer_type),
            !!data_col_name := str_remove_all(!!data_col_name, " |[:alpha:]|-|&|\\.$"),
            !!data_col_name := as.numeric(!!data_col_name)
        ) %>%
        select(cancer_type, !!data_col_name) %>%
        left_join(cancer_names, by = "cancer_type")
}


# SOURCE: https://cancerstatisticscenter.cancer.org/?_ga=2.9802422.205493334.1582111737-1270919539.1456557450#!/
ProjectTemplate::cache("acs_2020_estimates_df", {
    acs_2020_estimates_path <- file.path(
        DATA_DIR, "acs_NewCaseEstimates.xlsx"
    )

    acs_2020_estimates_df <- readxl::read_excel(acs_2020_estimates_path,
                                                sheet = "All US",
                                                skip = 6) %>%
        prepare_acs_data(est_incidence) %>%
        arrange(cancer, est_incidence)

    return(acs_2020_estimates_df)
})


# SOURCE: https://cancerstatisticscenter.cancer.org/?_ga=2.9802422.205493334.1582111737-1270919539.1456557450#!/
ProjectTemplate::cache("acs_incidence_df", {
    acs_incidence_path <- file.path(
        DATA_DIR, "acs_IncRate.xlsx"
    )

    acs_incidence_df <- readxl::read_excel(acs_incidence_path,
                                           sheet = "All US",
                                           skip = 6) %>%
        prepare_acs_data(incidence) %>%
        mutate(incidence = incidence * 1e5) %>%
        arrange(cancer, incidence)

    return(acs_incidence_df)
})
