# A second method for obtaining KRAS allele-specific synthetic lethal targets
# for Shikha.


GRAPHS_DIR <- "90_31_synlet-for-shikha_m2"
reset_graph_directory(GRAPHS_DIR)

TABLES_DIR <- GRAPHS_DIR
reset_table_directory(TABLES_DIR)

set.seed(0)

#### ---- Data preparation ---- ####

modeling_data <- depmap_modelling_df %>%
    filter_depmap_by_allele_count() %>%
    group_by(cancer) %>%
    filter(n_distinct(kras_allele) >= 3) %>%
    filter(!is_deleted) %>%
    add_count(cancer, hugo_symbol, kras_allele) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(all(n >= 3)) %>%
    select(-n) %>%
    ungroup()

modeling_data <- modeling_data %>%
    group_by(cancer, hugo_symbol) %>%
    mutate(rna_expression_std = scale_numeric(rna_expression,
                                              na.rm = TRUE)) %>%
    group_by(cancer, hugo_symbol) %>%
    nest() %>%
    ungroup()


#### ---- Model Experimentation ---- ####

set.seed(0)
eg_genes <- sample(unique(modeling_data$hugo_symbol), 20)
eg_genes
d <- modeling_data %>%
    filter(hugo_symbol %in% eg_genes)

d %>% count(cancer)
d %>%
    unnest(data) %>%
    write_tsv(table_path(TABLES_DIR, "sample-modeling-data.tsv"))
