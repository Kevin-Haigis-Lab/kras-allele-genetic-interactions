
info(logger, "Combining all RC-test results into a single data frame.")

cache("rc_test_results",
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

    read_rctest_output_file <- function(file_path) {
        stop("Good things come to those who STOP!")
    }

    all_files <- purrr::map_chr(
            c("comutation", "exclusivity"),
            ~ list.files(file.path("data", "rc-test", .x),
                         full.names = TRUE,
                         pattern = "rds")
        ) %>%
        unlist() %>%
        purrr::map(read_rctest_output_file) %>%
        bind_rows() %>%
        mutate(hugo_symbol = get_other_gene_from_genesets(gene_sets)) %>%
        filter(!is.na(hugo_symbol))

    return(rc_test_results)
})

