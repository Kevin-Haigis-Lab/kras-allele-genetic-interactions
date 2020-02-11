# Compare the association of genes in the genetic interactions with each other
# in the PPI.

GRAPHS_DIR <- "40_20_comut-dependency-genes-ppi-connectivity"
reset_graph_directory(GRAPHS_DIR)

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()


calc_nodes_connectivity <- function(from, to) {
    from_idx <- get_node_index(simple_combined_ppi_gr, name == !!from)
    to_idx <- get_node_index(simple_combined_ppi_gr, name == !!to)
    if (igraph::are_adjacent(simple_combined_ppi_gr, from, to)) {
        return (NA)
    }
    d <- igraph::vertex_connectivity(simple_combined_ppi_gr, from, to)
    return(as.numeric(d))
}
calc_nodes_connectivity <- memoise::memoise(calc_nodes_connectivity)



#### ---- Build null-distribution of connectivity ---- ####

ProjectTemplate::cache("null_connectivity_dist",
                       depends = "simple_combined_ppi_gr",
{
    set.seed(0)
    all_ppi_nodes <- na.omit(igraph::V(simple_combined_ppi_gr)$name)
    n_samples <- 1e3
    null_connectivity_dist <- tibble(
            from = sample(all_ppi_nodes, n_samples * 1.1, replace = TRUE),
            to = sample(all_ppi_nodes, n_samples * 1.1, replace = TRUE)
        ) %>%
            filter(from != to) %>%
            slice(1:n_samples) %>%
            mutate(connectivity = purrr::map2_dbl(from, to,
                                                  calc_nodes_connectivity)) %>%
            filter(is.finite(connectivity))
    return(null_connectivity_dist)
})

avg_null_connect <- mean(null_connectivity_dist$connectivity)
sd_null_connect <- sd(null_connectivity_dist$connectivity)

null_connectivity_hist <- null_connectivity_dist %>%
    ggplot() +
    geom_histogram(aes(x = connectivity), bins = 50,
                   fill = "grey75", color = "grey25", alpha = 0.5) +
    geom_vline(xintercept = avg_null_connect,
               linetype = 2, color = "dodgerblue", size = 0.6) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_text(hjust = 0.5)
    ) +
    labs(
        x = "shortest path length",
        y = "count",
        title = "Null distribution of shortest path lengths"
    )
ggsave_wrapper(
    null_connectivity_hist,
    plot_path(GRAPHS_DIR, "null_connectivity_hist.svg"),
    "small"
)


#### ---- Calculate  ---- ####


measure_connectivity <- function(df, status = NULL) {
    if (!is.null(status)) { print(status) }
    df %>%
        filter(hugo_symbol %in% igraph::V(simple_combined_ppi_gr)$name) %>%
        u_pull(hugo_symbol) %>%
        combn(2) %>%
        t() %>%
        as.data.frame() %>%
        as_tibble() %>%
        dplyr::rename(from = V1, to = V2) %>%
        mutate(
            connectivity = future_map2_dbl(from, to, calc_nodes_connectivity)
        )
}


library(furrr)
library(tictoc)
future::plan(multiprocess)

ProjectTemplate::cache("comut_dep_connectivity",
                       depends = c("depmap_gene_clusters_pairwise_df",
                                   "simple_combined_ppi_gr"),
{
    comut_dep_connectivity <- depmap_gene_clusters_pairwise_df %>%
        select(cancer, group1, group2) %>%
        unique() %>%
        group_by(cancer) %>%
        summarise(allele = list(unique(c(unlist(group1), unlist(group2))))) %>%
        ungroup() %>%
        unnest(allele) %>%
        filter(cancer == "COAD") %>%
        mutate(cancer_allele = paste(cancer, "-", allele)) %>%
        mutate(
            data = map2(cancer, allele, get_overlapped_df),
            connection_data = map2(data, cancer_allele, measure_connectivity)
        )
    return(comut_dep_connectivity)
})


comut_dep_connectivity %>%
    select(cancer, allele, connection_data) %>%
    unnest(connection_data) %>%
    filter(is.finite(connectivity)) %>%
    group_by(cancer, allele) %>%
    summarise(
        mean_connect = mean(connectivity),
        sd_connect = sd(connectivity)
    ) %>%
    ungroup()

comut_dep_connectivity_violins <- comut_dep_connectivity %>%
    select(cancer, allele, connection_data) %>%
    unnest(connection_data) %>%
    filter(is.finite(connectivity)) %>%
    mutate(allele = factor_alleles(allele)) %>%
    ggplot(aes(x = allele,  y = connectivity)) +
    facet_wrap(cancer ~ ., scales = "free") +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    scale_y_continuous(
        expand = expand_scale(mult = c(0, 0.02)),
        breaks = c(5, 10, 25, 50, 100, 150, 200, 250, 300)
    ) +
    coord_trans(y = my_trans_log10) +
    theme_bw(base_size = 7, base_family = "Arial")
ggsave_wrapper(
    comut_dep_connectivity_violins,
    plot_path(GRAPHS_DIR, "comut_dep_connectivity_violins.svg"),
    "medium"
)
