# Compare the association of genes in the genetic interactions with each other
# in the PPI.

GRAPHS_DIR <- "40_20_comut-dependency-genes-ppi-connectivity"
reset_graph_directory(GRAPHS_DIR)

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()



weighted_shortest_path <- function(gr, from, to, mode = "all") {
    paths <- igraph::all_shortest_paths(gr, from, to, mode = mode)
    n_paths <- length(paths$res)
    d <- igraph::distances(gr, from, to, mode = mode)
    return(d * sqrt(n_paths))
}

calc_nodes_connectivity <- function(from, to) {
    d <- igraph::distances(simple_combined_ppi_gr, from, to)
    return(as.numeric(d))
}
calc_nodes_connectivity <- memoise::memoise(calc_nodes_connectivity)



#### ---- Build null-distribution of connectivity ---- ####

ProjectTemplate::cache("null_connectivity_dist",
                       depends = "simple_combined_ppi_gr",
{
    set.seed(0)
    all_ppi_nodes <- na.omit(igraph::V(simple_combined_ppi_gr)$name)
    n_samples <- 2e3
    null_connectivity_dist <- tibble(
            from = sample(all_ppi_nodes, n_samples * 1.1, replace = TRUE),
            to = sample(all_ppi_nodes, n_samples * 1.1, replace = TRUE)
        ) %>%
            filter(from != to) %>%
            slice(1:n_samples) %>%
            mutate(connectivity = map2_dbl(from, to,
                                           calc_nodes_connectivity)) %>%
            filter(is.finite(connectivity))
    return(null_connectivity_dist)
})

avg_null_connect <- mean(null_connectivity_dist$connectivity)
sd_null_connect <- sd(null_connectivity_dist$connectivity)

null_connectivity_bar <- null_connectivity_dist %>%
    mutate(connectivity = as.character(connectivity)) %>%
    ggplot() +
    geom_bar(aes(x = '', fill = connectivity), position = "fill",
             alpha = 0.9) +
    scale_y_discrete(expand = expand_scale(mult = c(0, 0))) +
    scale_x_discrete(expand = expand_scale(mult = c(0, 0))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    labs(
        y = "geodesic distance",
        title = "Null dist. of\ngeodesic distance",
        fill = "distance"
    )
ggsave_wrapper(
    null_connectivity_bar,
    plot_path(GRAPHS_DIR, "null_connectivity_bar.svg"),
    width = 2, height = 4
)


#### ---- Calculate  ---- ####


measure_connectivity <- function(df, status = NULL, sample_at = 5e4) {
    if (!is.null(status)) { print(status) }

    res <- df %>%
        filter(hugo_symbol %in% igraph::V(simple_combined_ppi_gr)$name) %>%
        u_pull(hugo_symbol) %>%
        combn(2) %>%
        t() %>%
        as.data.frame() %>%
        as_tibble() %>%
        dplyr::rename(from = V1, to = V2)

    if (nrow(res) > sample_at) {
        cat(glue("(Sampling {sample_at} paths.)"), "\n")
        res <- sample_n(res, sample_at)
    }

    res %<>% mutate(connectivity = map2_dbl(from, to, calc_nodes_connectivity))

    return(res)
}

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

null_tib <- tibble(cancer = sort(unique(comut_dep_connectivity$cancer))) %>%
    mutate(
        connectivity = list(rep(null_connectivity_dist$connectivity, n()))
    ) %>%
    unnest(connectivity) %>%
    add_column(allele = "null dist.")

comut_dep_connectivity_bars <- comut_dep_connectivity %>%
    select(cancer, allele, connection_data) %>%
    unnest(connection_data) %>%
    filter(is.finite(connectivity)) %>%
    bind_rows(null_tib) %>%
    mutate(allele = factor(allele,
                           levels = c(names(short_allele_pal), "null dist.")),
           connectivity = as.character(connectivity)) %>%
    ggplot(aes(x = allele)) +
    facet_wrap(cancer ~ ., scales = "free", nrow = 1) +
    geom_bar(
        aes(fill = connectivity),
        position = "fill", alpha = 1.0, color = "black",
        width = 1.0
    ) +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 1) +
    scale_y_discrete(expand = expand_scale(mult = c(0, 0))) +
    scale_x_discrete(expand = expand_scale(mult = c(0, 0))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank()
    )

ggsave_wrapper(
    comut_dep_connectivity_bars,
    plot_path(GRAPHS_DIR, "comut_dep_connectivity_bars.svg"),
    "wide"
)

saveRDS(
    comut_dep_connectivity_bars,
    get_fig_proto_path("comut_dep_connectivity_bars", 5)
)
