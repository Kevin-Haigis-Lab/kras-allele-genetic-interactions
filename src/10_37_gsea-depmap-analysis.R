
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
    unnest(data) %>%
    mutate(
        gene_set_family = str_split_fixed(name, "_", 2)[, 1],
        gene_set = str_split_fixed(name, "_", 2)[, 2]
    )

cache("gsea_df")



#### ---- Plot GSEA results ---- ####


uninteresting_terms_regex <- c(
    "pancreas", "keratin", "disease", "muscle"
) %>%
    paste0(collapse = "|") %>%
    regex(ignore_case = TRUE)

plot_gsea_results <- function(cancer, gene_set_family, data) {
    mod_data <- data %>%
        filter(abs(nes) >= 1.2 & fdr_q_val < 0.2) %>%
        filter(!str_detect(gene_set, uninteresting_terms_regex))

    if (nrow(mod_data) == 0) { return() }

    p <- mod_data %>%
        mutate(gene_set = str_replace_us(gene_set),
               gene_set = str_to_sentence(gene_set)) %>%
        ggplot() +
        geom_point(
            aes(
                x = allele,
                y = gene_set,
                color = nes,
                size = -log10(fdr_q_val)
            )
        ) +
        scale_color_gradient2() +
        theme_bw() +
        theme(
            text = element_text("arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_blank(),
            axis.text.y = element_text(size = 10)
        ) +
        labs(
            title = glue("GSEA of Alleles in {cancer} ({gene_set_family})"),
            color = "NES",
            size = "-log10( adj. p-val. )"
        )
    save_path <- plot_path("10_37_gsea-depmap-analysis",
                           glue("gsea-results-{cancer}-{gene_set_family}.svg"))
    ggsave_wrapper(p, save_path, "wide")
}

gsea_df %>%
    filter(gene_set_family %in% c("HALLMARK", "KEGG", "REACTOME", "BIOCARTA", "PID")) %>%
    group_by(cancer, gene_set_family) %>%
    nest() %>%
    purrr::pwalk(plot_gsea_results)
