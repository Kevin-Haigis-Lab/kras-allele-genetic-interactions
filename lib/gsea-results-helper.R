# Functions for reading GSEA results.


# Standard filtering for GSEA results.
standard_gsea_results_filter <- function(df) {
  df %>%
    filter(abs(nes) >= 1.2 & fdr_q_val < 0.2) %>%
    filter(!str_detect(gene_set, uninteresting_terms_regex)) %>%
    mutate(fdr_q_val = purrr::map_dbl(fdr_q_val, ~ max(1 / 10000, .x)))
}

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
      unnest(data) %>%
      mutate(
        nes = -nes,
        es = -es
      )
    return(report_df)
  }
  purrr::map(dirs, get_gsea_report)
}


# Read in the GSEA results for a gene set.
# Pass the full file path to `xls_path`.
read_gsea_geneset_xls <- function(xls_path) {
  suppressWarnings(read_tsv(xls_path, col_types = cols())) %>%
    janitor::clean_names() %>%
    select(
      probe, rank_in_gene_list, rank_metric_score,
      running_es, core_enrichment
    ) %>%
    mutate(core_enrichment = core_enrichment == "Yes")
}
read_gsea_geneset_xls <- memoise::memoise(read_gsea_geneset_xls)


# Get the enrichment results for the gene set `name` in `cancer` and `allele`.
get_geneset_enrichment_results <- function(cancer, allele, name) {
  dir <- list.dirs(file.path("data", "gsea", "output"), recursive = FALSE)
  dir <- dir[str_detect(dir, cancer) & str_detect(dir, allele)]

  if (length(dir) != 1) {
    cat("Below are the directories:\n")
    print(dir)
    cat(glue("cancer: {cancer}, allele: {allele}"), "\n")
    stop("There is ambiguity about which directory to use.")
  }

  fpath <- list.files(dir, full.names = TRUE)
  idx <- basename(fpath) == paste0(name, ".xls")
  if (sum(idx) == 0) {
    all_file_names <- file_sans_ext(basename(fpath))
    idx <- purrr::map_lgl(all_file_names, ~ str_detect(.x, name))
  }
  fpath <- fpath[idx]

  if (length(fpath) > 1) {
    cat("Below are the file paths:\n")
    print(fpath)
    cat(glue("cancer: {cancer}, allele: {allele}"), "\n")
    cat("gene set:", name, "\n")
    stop("There is ambiguity about which file to use.")
  } else if (length(fpath) == 0) {
    cat("There are no files for the following parameters:\n")
    cat("    cancer:", cancer, "\n")
    cat("    allele:", allele, "\n")
    cat("  gene set:", name, "\n")
    return(NULL)
  }

  gsea_xls_data <- read_gsea_geneset_xls(fpath)

  # Reverse the ordering of the data frame if the first entry is not in the
  # core enrichment because that means the gene set was *negatively* enriched.
  if (!gsea_xls_data$core_enrichment[[1]]) {
    gsea_xls_data %<>% arrange(-rank_in_gene_list)
  }

  return(gsea_xls_data)
}
