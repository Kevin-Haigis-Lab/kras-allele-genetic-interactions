# Fundamental functions for Section 40 in "src" dealing with the integration
# of the comutation and genetic dependency results.


#### ---- Codify the clusters from synthetic lethal analysis ---- ####

# Get the average `gene_effect` for an allele.
# The function for calculate the average is in the argument `f`.
#   By default it is `mean()`.
mean_gene_effect <- function(allele, data, f = mean) {
    filter(data, allele == !!allele) %>% pull(gene_effect) %>% unlist() %>% f()
}


# Get the average `gene_effect` for all alleles except for `allele`.
# The function for calculate the average is in the argument `f`.
#   By default it is `mean()`.
mean_other_gene_effect <- function(allele, data, f = mean) {
    filter(data, allele != !!allele) %>% pull(gene_effect) %>% unlist() %>% f()
}


# Extract and parse the pair-wise comparisons from the linear modeling of
#   the synthetic lethal data.
extract_pw <- function(pw, df) {
    tidy_pw <- tidy(pw) %>%
        janitor::clean_names() %>%
        mutate(
            g1_avg = purrr::map_dbl(group1, mean_gene_effect, data = df),
            g2_avg = purrr::map_dbl(group2, mean_gene_effect, data = df),
            g1_other_avg = purrr::map_dbl(group1, mean_other_gene_effect,
                                          data = df),
            g2_other_avg = purrr::map_dbl(group2, mean_other_gene_effect,
                                          data = df)
        )
    return(tidy_pw)
}


# A tibble of the pairwise comparisons between all alleles for the genes in
#   the `depmap_gene_clusters`.
make_depmap_gene_clusters_pairwise_df <- function() {
    if (exists("depmap_gene_clusters_pairwise_df")) {
        cat("(`depmap_gene_clusters_pairwise_df` already exists)\n")
        return()
    }
    cat("Creating `depmap_gene_clusters_pairwise_df`.\n")

    depmap_gene_clusters_pairwise_df <- model1_tib %>%
        inner_join(depmap_gene_clusters, by = c("cancer", "hugo_symbol")) %>%
        mutate(
            allele_pairwise = purrr::map2(allele_pairwise, data, extract_pw)
        ) %>%
        select(cancer, hugo_symbol, gene_cls, allele_pairwise) %>%
        unnest(allele_pairwise) %>%
        dplyr::rename(adj_p_value = p_value)

    assign("depmap_gene_clusters_pairwise_df",
           depmap_gene_clusters_pairwise_df,
           envir = .GlobalEnv)
    ProjectTemplate::cache("depmap_gene_clusters_pairwise_df",
                           depends = c("model1_tib", "depmap_gene_clusters"))
    invisible()
}


prepare_simple_combined_ppi_gr <- function() {
    if (exists("simple_combined_ppi_gr")) {
        cat("(`simple_combined_ppi_gr` already exists)\n")
        return()
    }

    cat("Creating `simple_combined_ppi_gr`.\n")
    simple_combined_ppi_gr <- convert(combined_ppi_gr, to_simple) %E>%
        mutate(num_source = purrr::map_dbl(.orig_data,
                                           ~ n_distinct(.x$source))) %>%
        select(-.tidygraph_edge_index, -.orig_data) %N>%
        select(-.tidygraph_node_index)

    assign("simple_combined_ppi_gr",
           simple_combined_ppi_gr,
           envir = .GlobalEnv)
    ProjectTemplate::cache("simple_combined_ppi_gr",
                           depends = c("combined_ppi_gr"))
    invisible()
}



# Get the synthetic lethal data for a cancer and allele.
get_synthetic_lethal_data <- function(cancer, allele, adj_p_value = 0.05) {
    depmap_gene_clusters_pairwise_df %>%
        filter(cancer == !!cancer) %>%
        filter(adj_p_value < !!adj_p_value) %>%
        filter(group1 == !!allele | group2 == !!allele) %>%
        select(-cancer)
}


# Get the genetic interaction data for a cancer and allele.
get_genetic_interaction_data <- function(cancer, allele) {
    genetic_interaction_df %>%
        filter(cancer == !!cancer & allele == !!allele) %>%
        select(hugo_symbol, p_val, genetic_interaction) %>%
        mutate(genetic_interaction = ifelse(
            genetic_interaction == "exclusivity",
            "reduced\ncomut.", "increased\ncomut."
        ))
}


# Get a single data frame with all of the dependency and comutation
# data for a cancer and allele.
get_overlapped_df <- function(cancer, allele) {
    dependency_df <- get_synthetic_lethal_data(cancer, allele) %>%
        mutate(
            allele = !!allele,
            other_allele = ifelse(group1 == !!allele, group2, group1),
            comparison = paste(allele, "-", other_allele),
            allele_diff = ifelse(
                group1 == !!allele,
                g1_avg - g1_other_avg,
                g2_avg - g2_other_avg
            )
        ) %>%
        group_by(hugo_symbol) %>%
        slice(1) %>%
        ungroup()
    genetic_df <- get_genetic_interaction_data(cancer, allele)

    df <- full_join(dependency_df, genetic_df, by = "hugo_symbol") %>%
        mutate(interaction_source = case_when(
            is.na(gene_cls) ~ "comutation",
            is.na(genetic_interaction) ~ "dependency",
            TRUE ~ "both"
        ))
    return(df)
}


# Get a PPI annotated with comutation and genetic dependency data.
get_overlapped_gr <- function(cancer, allele, min_comp_size, ignore_genes) {
    merged_df <- get_overlapped_df(cancer, allele)
    gr <- simple_combined_ppi_gr %N>%
        filter(!(name %in% !!ignore_genes)) %>%
        filter(name %in% c(merged_df$hugo_symbol, "KRAS")) %>%
        jhcutils::filter_component_size(min_size = min_comp_size)

    return(list(graph = gr, data = merged_df))
}
