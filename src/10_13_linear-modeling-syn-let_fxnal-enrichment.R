
#### ---- Functional annotation of heatmap gene (row) clusters ---- ####

cluster_terms <- depmap_gene_clusters %>%
    group_by(cancer, gene_cls) %>%
    summarise(genes = list(hugo_symbol)) %>%
    ungroup() %>%
    mutate(enrichr_res = purrr::map(genes, enrichr_wrapper)) %>%
    select(-genes) %>%
    unnest(enrichr_res) %>%
    filter(!str_detect(term, !!uniteresting_enrichr_regex)) %>%
    mutate(n_genes = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.2 & n_genes >= 3)

cache("cluster_terms", depends = "model1_tib")


cluster_terms %>%
    filter(cancer == "PAAD") %>%
    select(datasource, term, adjusted_p_value, genes)

# Write the enrichment results to file.
cluster_terms %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "gene-clusters-fxnal-enrichment.tsv"
    ))

cluster_terms %>%
    group_by(cancer, term) %>%
    filter(n_distinct(gene_cls) == 1) %>%
    group_by(cancer, datasource, gene_cls) %>%
    arrange(adjusted_p_value, desc(n_genes)) %>%
    slice(1:5) %>%
    ungroup() %>%
    arrange(cancer, gene_cls, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "gene-clusters-fxnal-enrichment_top10.tsv"
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
    arrange(cancer, gene_cls, datasource, term, adjusted_p_value, desc(n_genes)) %>%
    write_tsv(
        file.path("tables",
                  "10_10_linear-modeling-syn-let",
                  "gene-clusters-fxnal-enrichment_uncommon.tsv"
    ))



