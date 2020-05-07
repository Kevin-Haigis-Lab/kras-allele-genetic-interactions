
GRAPHS_DIR <- "10_13_linear-modeling-syn-let_fxnal-enrichment"
reset_graph_directory(GRAPHS_DIR)

TABLES_DIR <- GRAPHS_DIR
reset_table_directory(TABLES_DIR)

#### ---- Functional annotation of heatmap gene (row) clusters ---- ####

cluster_terms <- depmap_gene_clusters %>%
    group_by(cancer, gene_cls) %>%
    summarise(genes = list(hugo_symbol)) %>%
    ungroup() %>%
    mutate(enrichr_res = purrr::map(genes, enrichr_wrapper)) %>%
    select(-genes) %>%
    unnest(enrichr_res) %>%
    filter(!str_detect(term, !!uninteresting_enrichr_regex)) %>%
    mutate(n_genes = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.2 & n_genes > 2) %>%
    mutate(datasource_hr = unlist(mapping_datasource_names[datasource]))

ProjectTemplate::cache("cluster_terms", depends = "depmap_gene_clusters")


# Write the enrichment results to file.
cluster_terms %>%
    write_tsv(table_path(TABLES_DIR, "gene-clusters-fxnal-enrichment.tsv"))

cluster_terms %>%
    group_by(cancer, term) %>%
    filter(n_distinct(gene_cls) == 1) %>%
    group_by(cancer, datasource, gene_cls) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:10) %>%
    ungroup() %>%
    arrange(cancer, gene_cls, datasource, term, adjusted_p_value,
            desc(n_genes)) %>%
    write_tsv(
        table_path(TABLES_DIR, "gene-clusters-fxnal-enrichment_top10.tsv"
    ))


common_term_regex <- c(
    "mitoch", "citric", "transcript", "translat"
) %>%
    paste0(collapse = "|") %>%
    regex(ignore_case = TRUE)

cluster_terms %>%
    filter(!str_detect(term, common_term_regex)) %>%
    group_by(cancer, term) %>%
    filter(n_distinct(gene_cls) == 1) %>%
    group_by(cancer, datasource, gene_cls) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:5) %>%
    ungroup() %>%
    arrange(cancer, gene_cls, datasource, term, adjusted_p_value,
            desc(n_genes)) %>%
    write_tsv(
        table_path(TABLES_DIR, "gene-clusters-fxnal-enrichment_uncommon.tsv"
    ))


#### ---- Functional enrichment for each cancer (not gene clusters) ---- ####

cluster_terms_cancer <- depmap_gene_clusters %>%
    group_by(cancer) %>%
    summarise(genes = list(hugo_symbol)) %>%
    ungroup() %>%
    mutate(enrichr_res = purrr::map(genes, enrichr_wrapper)) %>%
    select(-genes) %>%
    unnest(enrichr_res) %>%
    filter(!str_detect(term, !!uninteresting_enrichr_regex)) %>%
    mutate(n_genes = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.2 & n_genes > 2) %>%
    mutate(datasource_hr = unlist(mapping_datasource_names[datasource]))

cluster_terms_cancer %>%
    write_tsv(table_path(TABLES_DIR, "cancer-fxnal-enrichment.tsv"))

cluster_terms_cancer %>%
    group_by(cancer, term) %>%
    group_by(cancer, datasource) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:10) %>%
    ungroup() %>%
    arrange(cancer, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(table_path(TABLES_DIR, "cancer-fxnal-enrichment_top10.tsv"))

cluster_terms_cancer %>%
    filter(!str_detect(term, common_term_regex)) %>%
    group_by(cancer, datasource) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:10) %>%
    ungroup() %>%
    arrange(cancer, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(table_path(TABLES_DIR, "cancer-fxnal-enrichment_uncommon.tsv"))
