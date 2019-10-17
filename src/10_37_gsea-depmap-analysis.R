
# Analyzing the results of GSEA on the DepMap data.


#### ---- Read in results ---- ####

# Read in a GSEA report file.
#   They have the file extension "xls" but are really just TSV files.
#   `read_tsv()` adds a spare column that is removed before returning.
read_gsea_report_xls <- function(file_path) {
    df <- suppressWarnings(read_tsv(file_path, col_types = cols())) %>%
        janitor::clean_names() %>%
        select(-x12)
    return(df)
}


# Retrieve the GSEA report data.
get_gsea_reports <- function(dirs) {
    get_gsea_report <- function(dir) {
        report_df <- list.files(dir, pattern = "report.*xls", full.names = TRUE) %>%
            tibble(file_path = .) %>%
            mutate(
                gsea_group = str_split_fixed(basename(file_path), "_", 5)[, 4],
                data = map(file_path, ~ read_gsea_report_xls(.x))
            ) %>%
            unnest(data)
        return(report_df)
    }

    purrr::map(dirs, get_gsea_report)
}


gsea_df <- file.path("data", "gsea", "output") %>%
    list.dirs(recursive = FALSE) %>%
    tibble(dir = .) %>%
    mutate(
        dir_base = basename(dir),
        cancer = str_extract(dir_base, "^[:alpha:]+(?=_)"),
        allele = str_extract(dir_base, "(?<=_)[:alnum:]+(?=\\.G)"),
        timestamp = str_extract(dir_base, "(?<=Gsea\\.)[:digit:]+$"),
        timestamp = lubridate::as_datetime(as.numeric(timestamp) / 1000),
        data = get_gsea_reports(dir)
    ) %>%
    unnest(data)

cache("gsea_df")



