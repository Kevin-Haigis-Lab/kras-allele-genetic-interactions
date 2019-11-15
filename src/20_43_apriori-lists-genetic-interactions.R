
# Checking for genes-of-interest (goi) in the genetic interactions with KRAS

library(ggraph)

#### ---- Manual inspection of hits ---- ####

cache("wide_genetic_interaction_df", depends = "genetic_interaction_gr",
{
    wide_genetic_interaction_df <- kegg_geneset_df %>%
        group_by(hugo_symbol) %>%
        summarise(KEGG = paste(gene_set, collapse = ", ")) %>%
        ungroup() %>%
        full_join(
            {
                cosmic_cgc_df %>%
                    select(hugo_symbol) %>%
                    add_column(CGC = TRUE)
            },
            by = "hugo_symbol"
        ) %>%
        full_join(
            {
                kras_interactors_bioid_df %>%
                    select(hugo_symbol) %>%
                    add_column(BioID = TRUE)
            },
            by = "hugo_symbol"
        )
    return(wide_genetic_interaction_df)
})



plot_genetic_interaction_graph <- function(gr_to_plot, CANCER, SUFFIX = "") {
    set.seed(0)
    gr_plot <- gr_to_plot %>%
        ggraph(layout = "stress") +
        geom_edge_link(aes(color = genetic_interaction, width = -log(p_val + 0.0000001))) +
        scale_edge_color_manual(values = comut_mutex_pal) +
        scale_edge_width_continuous(range = c(0.2, 1.5)) +
        geom_node_point(aes(color = node_color, size = node_size)) +
        scale_color_manual(values = short_allele_pal, na.value = "grey75") +
        scale_size_manual(values = c(big = 2, small = 1), guide = FALSE) +
        geom_node_text(aes(label = node_label), repel = TRUE, family = "Arial", size = 2) +
        theme_graph() +
        theme(
            text = element_text(family = "Arial")
        ) +
        labs(
            title = glue("Genes of interest with genetic interactions\nwith KRAS in {CANCER}"),
            color = "KRAS allele",
            edge_color = "genetic\ninteraction",
            edge_width = "-log( p-value )"
        )
    save_path <- plot_path(
        "20_43_apriori-lists-genetic-interactions",
        glue("goi_overlap_genetic_interactions_network_{CANCER}{SUFFIX}.svg")
    )
    ggsave_wrapper(gr_plot, save_path, "wide")
}



for (CANCER in unique(genetic_interaction_df$cancer)) {
    gr_to_plot <- genetic_interaction_gr %N>%
        left_join(wide_genetic_interaction_df, by = c("name" = "hugo_symbol")) %E>%
        filter(cancer == !!CANCER) %N>%
        filter(is_kras | !is.na(KEGG) | !is.na(CGC) | !is.na(BioID)) %>%
        filter(centrality_degree(mode = "all") > 0) %>%
        mutate(node_label = str_remove_all(name, "KRAS_"),
               node_color = ifelse(is_kras, node_label, NA),
               node_size = ifelse(is_kras, "big", "small"))

    if (igraph::vcount(gr_to_plot) == 0) { next }

    # plot all interactions with goi
    plot_genetic_interaction_graph(gr_to_plot, CANCER, "_allLists")

    # only KEGG
    gr_to_plot_MOD <- gr_to_plot %N>%
        filter(is_kras | !is.na(KEGG)) %>%
        filter(centrality_degree(mode = "all") > 0)
    if (igraph::vcount(gr_to_plot_MOD) > 0) {
        plot_genetic_interaction_graph(gr_to_plot_MOD, CANCER, "_kegg")
    }

    # only CGC
    gr_to_plot_MOD <- gr_to_plot %N>%
        filter(is_kras | !is.na(CGC)) %>%
        filter(centrality_degree(mode = "all") > 0)
    if (igraph::vcount(gr_to_plot_MOD) > 0) {
        plot_genetic_interaction_graph(gr_to_plot_MOD, CANCER, "_cgc")
    }

    # only BioID
    gr_to_plot_MOD <- gr_to_plot %N>%
        filter(is_kras | !is.na(BioID)) %>%
        filter(centrality_degree(mode = "all") > 0)
    if (igraph::vcount(gr_to_plot_MOD) > 0) {
        plot_genetic_interaction_graph(gr_to_plot_MOD, CANCER, "_BioID")
    }

}
