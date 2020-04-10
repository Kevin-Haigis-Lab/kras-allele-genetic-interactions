# Are the comutation genes for each KRAS allele closer to each other in the PPI?

GRAPHS_DIR <- "90_13_comutation-network_withingroup-ppi-closeness"
reset_graph_directory(GRAPHS_DIR)



# Measure the closeness of two genes, `from` and `to`, in a graph `gr`.
# `from` and `to` should be node names.
calculate_closeness <- function(gr, from, to) {
    from_i <- jhcutils::get_node_index(gr, name == from)
    to_i <- jhcutils::get_node_index(gr, name == to)
    d <- igraph::distances(gr, from_i, to_i)
    return(d)
}
# calculate_closeness <- memoise::memoise(calculate_closeness)


# Returns a tibble of all the possible combinations of genes in `gs`.
# The columns of the output tibble are "from" and "to".
combination_tibble <- function(gs, m = 2) {
    combn(gs, m, stringsAsFactors = FALSE) %>%
        t() %>%
        as.data.frame() %>%
        as_tibble() %>%
        dplyr::rename(from = V1, to = V2)
}


# Measures the closeness of all the genes in `gs` to each other in the PPI `gr`.
measure_ppi_closeness <- function(gs, gr, max_val = Inf) {
    gs <- gs[gs %in% igraph::V(gr)$name]
    if (length(gs) < 2) { return(NULL) }

    gene_combs <- combination_tibble(gs) %>%
        mutate(
            distance = map2_dbl(from, to, calculate_closeness, gr = gr),
            distance = ifelse(is.infinite(distance), max_val, distance)
        )
    return(gene_combs)
}


# Returns a list of the genes in the comutation network of an allele for
# a cancer.
get_genetic_interactions <- function(allele, cancer) {
    genetic_interaction_df %>%
        filter(cancer == !!cancer & allele %in% !!allele) %>%
        pull(hugo_symbol) %>%
        unique()
}
get_genetic_interactions <- memoise::memoise(get_genetic_interactions)


DIAMETER_OF_GIANT_COMP <- string_gr %>%
        jhcutils::get_giant_component() %>%
        igraph::diameter(directed = FALSE)


ProjectTemplate::cache("comutation_closeness_df",
                       depends = "genetic_interaction_df",
{
    comutation_closeness_df <- genetic_interaction_df %>%
        select(cancer, allele) %>%
        dplyr::rename(allele1 = allele) %>%
        mutate(allele2 = allele1) %>%
        group_by(cancer) %>%
        tidyr::complete(allele1, allele2) %>%
        ungroup() %>%
        unique() %>%
        mutate(
            genes1 = purrr::map2(allele1, cancer, get_genetic_interactions),
            genes2 = purrr::map2(allele2, cancer, get_genetic_interactions),
            genes = purrr::map2(genes1, genes2, ~ unique(unlist(c(.x, .y))))
        ) %>%
        select(-genes1, -genes2) %>%
        filter(cancer == "COAD") %>%
        mutate(
            closeness_scores = purrr::map(genes, measure_ppi_closeness,
                                          gr = string_gr,
                                          max_val = DIAMETER_OF_GIANT_COMP + 1)
        )
    return(comutation_closeness_df)
})



#### ---- Results analysis ---- ####

comutation_closeness_df %>%
    unnest(closeness_scores) %>%
    group_by(cancer, allele1, allele2) %>%
    summarise(
        avg_closeness = mean(distance),
        sd_closeness = sd(distance)
    ) %>%
    ungroup()


closeness_violins <- comutation_closeness_df %>%
    filter(cancer == "COAD") %>%
    unnest(closeness_scores) %>%
    ggplot(aes(x = allele2, y = distance)) +
    facet_grid(allele1 ~ .) +
    geom_boxplot(aes(fill = allele2), alpha = 0.5) +
    scale_fill_manual(values = short_allele_pal, guide = FALSE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank()
    ) +
    labs(
        x = "KRAS allele 1",
        y = "distances between genes",
        title = "Distance between genes in the COAD comutation network")
ggsave_wrapper(
    closeness_violins,
    plot_path(GRAPHS_DIR, "closeness_violins_COAD.svg"),
    width = 6, height = 12
)


closeness_heatmap <- comutation_closeness_df %>%
    filter(cancer == "COAD") %>%
    unnest(closeness_scores) %>%
    group_by(cancer, allele1, allele2) %>%
    summarise(
        avg_closeness = mean(distance),
        sd_closeness = sd(distance)
    ) %>%
    group_by(cancer, allele1) %>%
    mutate(scaled_avg_closeness = scale(avg_closeness)[, 1]) %>%
    ggplot(aes(x = allele1, y = allele2)) +
    geom_tile(aes(fill = scaled_avg_closeness)) +
    scale_fill_gradient2() +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "Allele", y = "Comparison")
ggsave_wrapper(
    closeness_heatmap,
    plot_path(GRAPHS_DIR, "closeness_heatmap_COAD.svg"),
    "medium"
)


## CONCLUSION: see notebook from 2020-02-04
# https://project-notebook-comutation.netlify.com/post/2020-02-04/
