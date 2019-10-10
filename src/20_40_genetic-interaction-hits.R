
library(ggraph)


gene_tib <- genes_of_interest_df %>%
    group_by(hugo_symbol) %>%
    summarise(sources = paste(source, collapse = ", "))

gi_goi_overlap <- genetic_interaction_df %>%
    left_join(gene_tib, by = "hugo_symbol") %>%
    filter(!is.na(sources))

coad_gi_goi_overlap <- gi_goi_overlap %>% filter(cancer == "COAD" & allele == "G13D")

genes <- coad_gi_goi_overlap$hugo_symbol %>%
    unlist() %>%
    c("KRAS") %>%
    unique()
length(genes)

gr_plot <- string_gr %N>%
    filter(name %in% !!genes) %E>%
    filter(combined_score > 700) %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    ggraph(layout = "stress") +
    geom_edge_link(aes(width = combined_score, alpha = combined_score), color = "grey50") +
    scale_edge_width_continuous(range = c(0.1, 0.7)) +
    geom_node_point(size = 2, color = "grey75") +
    geom_node_text(aes(label = name), size = 2.5, repel = TRUE, family = "Arial", color = "black") +
    theme_graph() +
    theme(
        text = element_text(family = "arial")
    )
save_path <- plot_path("20_40_genetic-interaction-hits", "TEST_gr_coad.svg")
ggsave_wrapper(gr_plot, save_path, "large")



for (CANCER in unique(gi_goi_overlap$cancer)) {
    ALLELES <- gi_goi_overlap %>% filter(cancer == !!CANCER) %>% pull(allele) %>% unique()
    for (ALLELE in ALLELES) {
        genes <- gi_goi_overlap %>%
            filter(cancer == !!CANCER & allele == !!ALLELE) %>%
            pull(hugo_symbol) %>%
            unlist() %>%
            c("KRAS") %>%
            unique()

        gr <- string_gr %N>%
            mutate(
                keeper = name %in% !!genes,
                neighbor_of_keeper = map_local_int(
                    .f = num_qual_neighbors,
                    lgl_filter = rlang::expr(keeper)
                )
            ) %>%
            filter(keeper | neighbor_of_keeper)
    }
}