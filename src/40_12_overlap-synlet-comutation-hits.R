# Individual genes that have both comutation and genetic dpendency interactions
# with an allele in a cancer.

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()

#### ---- Overlap of comutation and genetic dependency analysis ---- ####

get_shared_comutation_dependency <- function(cancer, allele) {
    get_overlapped_df(cancer, allele) %>%
        select(-c(group1, group2, g1_avg, g2_avg,
                  g1_other_avg, g2_other_avg)) %>%
        filter(interaction_source == "both") %>%
        mutate(cancer = !!cancer) %>%
        select(cancer, allele, everything())
}


comut_dep_overlap_tbl <- depmap_gene_clusters_pairwise_df %>%
    select(cancer, group1, group2) %>%
    unique() %>%
    group_by(cancer) %>%
    summarise(allele = list(unique(c(unlist(group1), unlist(group2))))) %>%
    ungroup() %>%
    unnest(allele) %>%
    pmap(get_shared_comutation_dependency) %>%
    bind_rows() %>%
    mutate(genetic_interaction = str_replace(genetic_interaction,
                                             "\\ncomutation", " comut.")) %>%
    select(cancer, allele, hugo_symbol,
           comparison, adj_p_value,
           genetic_interaction, p_val)

knitr::kable(comut_dep_overlap_tbl)

saveRDS(
    comut_dep_overlap_tbl,
    get_fig_proto_path("comut_dep_overlap_tbl.rds", 5)
)
