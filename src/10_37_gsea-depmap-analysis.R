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


# Pull a waterfall plot ("enplot"?) from the GSEA results for a the specified
#   cancer, allele, and name (gene-set).
pull_gsea_enplot <- function(cancer, allele, gene_set, gene_set_family,
                             force_overwrite = FALSE,
                             ...) {
    cancer_allele <- paste0(cancer, "_", allele)

    gsea_dirs <- file.path("data", "gsea", "output") %>%
        list.dirs(recursive = FALSE)
    gsea_dir <- gsea_dirs[str_detect(gsea_dirs, cancer_allele)]

    svg_files <- list.files(gsea_dir, pattern = "^enplot.*svg.gz", full.names = TRUE)

    if (str_length(gene_set) > 100) {
        gene_set_regex <- glue("{gene_set_family}_{str_sub(gene_set, 1, 100)}")
    } else {
        gene_set_regex <- glue("{gene_set_family}_{gene_set}_[:digit:]+")
    }

    svg_file <- svg_files[str_detect(basename(svg_files), gene_set_regex)]

    target_dir <- file.path("graphs", "10_37_gsea-depmap-output", cancer_allele)
    if (!dir.exists(target_dir)) {
        dir.create(target_dir)
    }

    target_path <- file.path(target_dir, basename(svg_file))
    unziped_target_path <- gsub("[.]gz$", "", target_path)

    if (length(target_path) != 1) { return(NULL) }

    if (file.exists(unziped_target_path) & !force_overwrite) {
        return(NULL)
    }

    a <- file.copy(svg_file, target_path)
    a <- R.utils::gunzip(target_path, overwrite = TRUE)

    return(NULL)
}


#### ---- Plot GSEA results ---- ####


uninteresting_terms_regex <- c(
    "pancreas", "keratin", "disease", "muscle"
) %>%
    paste0(collapse = "|") %>%
    regex(ignore_case = TRUE)


gsea_plot <- function(tib, title_suffix = "") {
    p <- tib %>%
        mutate(gene_set = str_replace_us(gene_set),
               gene_set = str_to_sentence(gene_set),
               gene_set = str_wrap(gene_set, width = 50)) %>%
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
            axis.text.y = element_text(size = 7)
        ) +
        labs(
            title = glue("GSEA of Alleles in {title_suffix}"),
            color = "NES",
            size = "-log10( adj. p-val. )"
        )
    return(p)
}


plot_gsea_results <- function(cancer, data) {
    mod_data <- data %>%
        filter(abs(nes) >= 1.2 & fdr_q_val < 0.2) %>%
        filter(!str_detect(gene_set, uninteresting_terms_regex))

    if (nrow(mod_data) == 0) { return() }

    purrr::pwalk(mod_data, pull_gsea_enplot, cancer = cancer)

    p <- gsea_plot(mod_data, title_suffix = cancer)
    save_path <- plot_path("10_37_gsea-depmap-analysis",
                           glue("gsea-results-{cancer}-all.svg"))
    ggsave_wrapper(p, save_path, "large")

    for (fam in unique(data$gene_set_family)) {
        suffix <- glue("{cancer} ({fam})")
        tt <- mod_data %>%
            filter(gene_set_family == !!fam)

        if (nrow(tt) == 0) { next }

        p <- gsea_plot(tt, title_suffix = suffix)
        save_path <- plot_path("10_37_gsea-depmap-analysis",
                               glue("gsea-results-{cancer}-{fam}.svg"))
        ggsave_wrapper(p, save_path, "large")
    }
}

gsea_df %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    filter(gene_set_family %in% c("HALLMARK", "KEGG", "REACTOME", "BIOCARTA", "PID")) %>%
    group_by(cancer) %>%
    nest() %>%
    purrr::pwalk(plot_gsea_results)



#### ---- Plot select GSEA results ---- ####


select_gsea_results <- list(
    COAD = c(
        "vegfr2 mediated cell proliferation",
        "tp53 regulates metabolic genes",
        "srp dependent cotranslational protein targeting to membrane",
        "respiratory electron transport",
        "regulation of expression of slits and robos",
        "complement cascade",
        "oxidative phosphorylation",
        "gpcr pathway",
        "nonsense mediated decay nmd"
    ),
    LUAD = c(
        "srp dependent cotranslational protein targeting to membrane",
        "nkt pathway",
        "fanconi anemia pathway",
        "eukaryotic translation initiation",
        "hdr through homologous recombination hrr",
        "beta alanine metabolism",
        "bard1 pathway",
        "nonsense mediated decay nmd",
        "steroid hormone biosynthesis"
    ),
    PAAD = c(
        "toll pathway",
        "tp53 regulates metabolic genes",
        "regulation of cholesterol biosynthesis by srebp srebf",
        "tgfb pathway",
        "nfkb pathway",
        "jnk c jun kinases phosphorylation and activation mediated by activated human tak1",
        "hedgehog signaling",
        "g2 m dna damage checkpoint",
        "fak pathway"
    )
)


standardize_names <- function(x) {
    str_to_lower(x) %>%
        str_replace_all("_", " ")
}

filter_gsea_for_select_results <- function(cancer, df, key) {
    gene_sets_to_keep <- standardize_names(unlist(key[[cancer]]))
    df %>%
        mutate(.gene_set = standardize_names(gene_set)) %>%
        filter(.gene_set %in% !!gene_sets_to_keep) %>%
        select(-.gene_set)
}


select_gsea_plot <- function(cancer, data, ...) {
    mod_data <- data %>%
        filter(abs(nes) >= 1.2 & fdr_q_val < 0.2) %>%
        filter(!str_detect(gene_set, uninteresting_terms_regex))

    if (nrow(mod_data) == 0) { return() }

    p <- gsea_plot(mod_data, title_suffix = cancer)
    save_path <- plot_path("10_37_gsea-depmap-analysis",
                           glue("gsea-results-{cancer}-select.svg"))
    ggsave_wrapper(p, save_path, "wide")
}


gsea_df %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    filter(gene_set_family %in% c("HALLMARK", "KEGG", "REACTOME", "BIOCARTA", "PID")) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(data = purrr::map2(cancer, data,
                              filter_gsea_for_select_results,
                              key = select_gsea_results)) %>%
    purrr::pwalk(select_gsea_plot)