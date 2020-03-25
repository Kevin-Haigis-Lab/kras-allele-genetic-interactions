
# Checking for genes-of-interest (goi) in the genetic interactions with KRAS

GRAPHS_DIR <- "20_43_apriori-lists-genetic-interactions"
reset_graph_directory(GRAPHS_DIR)


#### ---- Manual inspection of hits ---- ####

ProjectTemplate::cache("wide_genetic_interaction_df",
                       depends = "genetic_interaction_gr",
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


# Queryable list of which plots to save to protos for Figure 2.
imgs_to_save_for_figure <- list(
    COAD = list(suffix = "_allLists", num = 2, supp = FALSE),
    PAAD = list(suffix = "_allLists", num = 3, supp = FALSE)
)


plot_genetic_interaction_graph <- function(gr_to_plot, CANCER, SUFFIX = "") {
    set.seed(0)
    gr_plot <- gr_to_plot %>%
        ggraph(layout = "nicely") +
        geom_edge_link(
            aes(color = genetic_interaction,
                width = -log(p_val + 0.0000001))
        ) +
        scale_edge_color_manual(
            values = comut_mutex_pal,
            guide = FALSE
        ) +
        scale_edge_width_continuous(
            range = c(0.2, 1.5),
            guide = guide_legend(
                label.position = "top",
                keyheight = unit(1, "mm")
            )
        ) +
        geom_node_point(
            aes(color = node_color,
                size = node_size)
        ) +
        scale_color_manual(
            values = short_allele_pal,
            guide = FALSE,
            na.value = "grey75"
        ) +
        scale_size_manual(
            values = c(big = 2, small = 1),
            guide = FALSE
        ) +
        geom_node_text(
            aes(label = node_label,
                fontface = label_face),
            repel = TRUE,
            family = "Arial",
            size = 2
        ) +
        theme_graph() +
        theme(
            text = element_text(family = "Arial")
        ) +
        labs(
            edge_width = "-log( p-value )"
        )
    save_path <- plot_path(
        GRAPHS_DIR,
        glue("goi_overlap_genetic_interactions_network_{CANCER}{SUFFIX}.svg")
    )
    ggsave_wrapper(gr_plot, save_path, "wide")

    fig_info <- imgs_to_save_for_figure[[CANCER]]
    if (!is.null(fig_info) & SUFFIX %in% fig_info$suffix) {
        base_n <- file_sans_ext(basename(save_path))
        saveFigRds(gr_plot, base_n)
    }

}



for (CANCER in unique(genetic_interaction_df$cancer)) {
    gr_to_plot <- genetic_interaction_gr %N>%
        left_join(wide_genetic_interaction_df,
                  by = c("name" = "hugo_symbol")) %E>%
        filter(cancer == !!CANCER) %N>%
        filter(is_kras | !is.na(KEGG) | !is.na(CGC) | !is.na(BioID)) %>%
        filter(centrality_degree(mode = "all") > 0) %>%
        mutate(
            label_face = ifelse(is_kras, "bold", "plain"),
            node_label = str_remove_all(name, "KRAS_"),
            node_color = ifelse(is_kras, node_label, NA),
            node_size = ifelse(is_kras, "big", "small")
        )

    if (igraph::vcount(gr_to_plot) == 0) { next }

    # plot all interactions with goi
    plot_genetic_interaction_graph(gr_to_plot, CANCER, "_allLists")

    plot_specific_lists <- function(gene_list, suffix) {
        j <- which(colnames(as_tibble(gr_to_plot, "nodes")) == gene_list)
        idx <- !is.na(as_tibble(gr_to_plot, "nodes")[, j])

        gr_to_plot_MOD <- gr_to_plot %N>%
            filter(is_kras | !!idx) %>%
            filter(centrality_degree(mode = "all") > 0)
        if (igraph::vcount(gr_to_plot_MOD) > 0) {
            plot_genetic_interaction_graph(gr_to_plot_MOD, CANCER, suffix)
        }
    }

    tibble(gene_list = c("KEGG", "CGC", "BioID"),
           suffix = c("_kegg", "_cgc", "_BioID")) %>%
           pmap(plot_specific_lists)
}
