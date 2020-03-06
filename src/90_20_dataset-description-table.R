# Create a table with information about the data sources.

TABLES_DIR <- "90_20_dataset-description-table"
reset_table_directory(TABLES_DIR)

# Turn an array/list into a single string with comma separations.
list_to_csv_string <- function(a) {
    paste0(a, collapse = ", ")
}


#### ---- Tables of data source summaries ---- ####

# Count the number of unique hypermutant samples.
count_hypermuts <- function(tsb, is_hm) {
    n_distinct(tsb[is_hm])
}


# Data frame with summary statistics for each data source.
datasource_summary_df <- cancer_full_muts_df %>%
    filter(cancer != "SKCM") %>%
    group_by(cancer, dataset, target) %>%
    summarize(
        num_cancer_samples = n_distinct(tumor_sample_barcode),
        unique_genes_detected = n_distinct(hugo_symbol),
        types_of_mutations = list_to_csv_string(unique(mutation_type)),
        num_hypermutants = count_hypermuts(tumor_sample_barcode, is_hypermutant)
    ) %>%
    ungroup() %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    arrange(cancer, -num_cancer_samples) %T>%
    write_tsv(table_path(TABLES_DIR, "datasource-summaries.tsv"))


print_pretty_samplecount_tibble <- function(df, with_marginal = TRUE) {
    col_names <- c("cancer", "total samples", "total samples (no hypermutants)")
    df %>%
        group_by(cancer) %>%
        summarise(
            total_samples = sum(num_cancer_samples),
            total_samples_no_hypermut = total_samples - sum(num_hypermutants)
        ) %>%
        janitor::adorn_totals("row") %>%
        knitr::kable(col.names = col_names)
}

# Print out total number of samples per cancer.
cat("Total number of samples per cancer.\n")
print_pretty_samplecount_tibble(datasource_summary_df)

# Print out number of samples for exome/genome.
cat("Number of WXS or WGS samples per cancer.\n")
datasource_summary_df %>%
    filter(target %in% c("exome", "genome")) %>%
    print_pretty_samplecount_tibble()

# Print out number of samples for panel.
cat("Number of targeted sequencing samples per cancer.\n")
datasource_summary_df %>%
    filter(!target %in% c("exome", "genome")) %>%
    print_pretty_samplecount_tibble()



#### ---- Targeted panel size range ---- ####

# Data frame with summary of panel data sources.
paneldata_summary_df <- cancer_full_muts_df %>%
    filter(cancer != "SKCM") %>%
    filter(!target %in% c("genome", "exome")) %>%
    group_by(cancer, dataset, target) %>%
    summarize(
        num_cancer_samples = n_distinct(tumor_sample_barcode),
        num_genes_detected = n_distinct(hugo_symbol),
        genes_detected = list_to_csv_string(sort(unique(hugo_symbol)))
    ) %>%
    ungroup() %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    arrange(cancer, -num_cancer_samples, -num_genes_detected, dataset) %T>%
    write_tsv(table_path(TABLES_DIR, "panel-datasource-summaries.tsv"))

paneldata_summary_df %>%
    pull(num_genes_detected) %>%
    summary()
