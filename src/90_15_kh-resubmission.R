
# Code for generating the information and figures for KH's R01 resubmission.

library(ggraph)

# For reproducibility.
set.seed(0)

# Set-up directory to which images are saved.
SAVE_DIR <- "90_15_kh-resubmission"
reset_graph_directory(SAVE_DIR)

# Parameterize the cancer in case there are more in the future.
CANCER <- "COAD"

# Genes of interest.
goi_genes <- unique(genes_of_interest_df$hugo_symbol)

make_thick_network_plot <- function(gr) {
    p <- gr %>%
        ggraph(layout = "nicely") +
        geom_edge_link(
            aes(color = genetic_interaction),
            width = 0.3
        ) +
        scale_edge_color_manual(values = comut_mutex_pal) +
        geom_node_point(
            aes(color = node_color,
                size = node_size),
        ) +
        geom_node_text(
            aes(label = node_label,
                size = label_size),
            repel = TRUE,
            family = "Arial"
        ) +
        scale_color_manual(values = short_allele_pal, na.value = NA) +
        scale_size_identity() +
        theme_graph(base_family = "Arial", base_size = 8) +
        theme(
            plot.title = element_blank(),
            legend.position = "none"
        )
    return(p)
}



gr_plot <- genetic_interaction_gr %E>%
    filter(cancer == !!CANCER) %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(node_label = ifelse(is_kras, name, NA),
           node_label = str_remove_all(node_label, "KRAS_"),
           node_color = node_label,
           node_size = 2,
           label_size = 5) %>%
    make_thick_network_plot() +
    labs(
        title = glue("Genetic interactions in {CANCER}"),
        color = "KRAS allele",
        edge_color = "genetic\ninteraction"
    )
save_path <- plot_path(SAVE_DIR,
                       glue("genetic_interaction_network_{CANCER}_thick.svg"))
ggsave_wrapper(gr_plot, save_path, size = "medium")


gr_plot <- genetic_interaction_gr %E>%
    filter(cancer == !!CANCER & genetic_interaction == "comutation") %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(
        node_label = ifelse(is_kras | name %in% !!goi_genes, name, NA),
        node_label = str_remove_all(node_label, "KRAS_"),
        node_color = case_when(
            is_kras ~ node_label,
            name %in% !!goi_genes ~ "Other"
        ),
        node_size = case_when(
            is_kras ~ 2,
            name %in% !!goi_genes ~ 0.5
        ),
        label_size = case_when(
            is_kras ~ 5,
            name %in% !!goi_genes ~ 2.5
        )
    ) %>%
    make_thick_network_plot() +
    labs(
        title = glue("Genetic interactions in {CANCER}"),
        color = "KRAS allele",
        edge_color = "genetic\ninteraction"
    )

save_path <- plot_path(SAVE_DIR,
                       glue("genetic_interaction_network_{CANCER}_thick_comutation_goi.svg"))
ggsave_wrapper(gr_plot, save_path, size = "medium")


