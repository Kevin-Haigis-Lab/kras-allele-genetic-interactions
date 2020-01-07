# Create a table with information about the data sources.

TABLES_DIR <- "90_20_dataset-description-table"
reset_table_directory(TABLES_DIR)

# Turn an array/list into a single string with comma separations.
list_to_csv_string <- function(a) {
    paste0(a, collapse = ", ")
}


#### ---- Tables of data source summaries ---- ####


# Data frame with summary statistics for each data source.
datasource_summary_df <- cancer_full_muts_df %>%
    filter(cancer != "SKCM") %>%
    group_by(cancer, dataset, target) %>%
    summarize(
        num_cancer_samples = n_distinct(tumor_sample_barcode),
        unique_genes_detected = n_distinct(hugo_symbol),
        types_of_mutations = list_to_csv_string(unique(mutation_type))
    ) %>%
    ungroup() %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    arrange(cancer, -num_cancer_samples) %T>%
    write_tsv(table_path(TABLES_DIR, "datasource-summaries.tsv"))


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


#### ---- GENIE panel size range ---- ####

paneldata_summary_df %>%
    pull(num_genes_detected) %>%
    summary()
