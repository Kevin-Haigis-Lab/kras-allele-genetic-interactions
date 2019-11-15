# Make a data structure for the overlap of the synthetic lethal and genetic
# interaction analyses.

library(ggraph)

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
            g1_other_avg = purrr::map_dbl(group1, mean_other_gene_effect, data = df),
            g2_other_avg = purrr::map_dbl(group2, mean_other_gene_effect, data = df)
        )
    return(tidy_pw)
}


# A tibble of the pairwise comparisons between all alleles for the genes in
#   the `depmap_gene_clusters`.
depmap_gene_clusters_pairwise_df <- model1_tib %>%
    inner_join(depmap_gene_clusters, by = c("cancer", "hugo_symbol")) %>%
    mutate(allele_pairwise = purrr::map2(allele_pairwise, data, extract_pw)) %>%
    select(cancer, hugo_symbol, gene_cls, allele_pairwise) %>%
    unnest(allele_pairwise) %>%
    dplyr::rename(adj_p_value = p_value)


# Get the synthetic lethal data for a cancer and allele.
get_synthetic_lethal_data <- function(cancer, allele) {
    depmap_gene_clusters_pairwise_df %>%
        filter(cancer == !!cancer) %>%
        filter(adj_p_value < 0.10) %>%
        filter(group1 == !!allele | group2 == !!allele) %>%
        select(-cancer)
}

# Get the genetic interaction data for a cancer and allele.
get_genetic_interaction_data <- function(cancer, allele) {
    genetic_interaction_df %>%
        filter(cancer == !!cancer & allele == !!allele) %>%
        select(hugo_symbol, p_val, genetic_interaction)
}

# A filtered STRING PPI network
string_gr_prefilt <- string_gr %E>% filter(combined_score >= 800)


network_node_shapes <- list(
    "genetic" = 15,
    "synthetic lethal" = 17,
    "both" = 18,
    "other" = 19
)

network_node_colors <- list(
    "comutation" = comut_mutex_pal[["comutation"]],
    "exclusivity" = comut_mutex_pal[["exclusivity"]],
    "synthetic lethal up" = synthetic_lethal_pal[["up"]],
    "synthetic lethal down" = synthetic_lethal_pal[["down"]],
    "both" = "#ED6165",
    "other" = "#F44C4C"
) %>%
    lighten(factor = 1.4)

# Make the purple even brighter.
network_node_colors["synthetic lethal down"] <- lighten(
    network_node_colors[["synthetic lethal down"]],
    factor = 1.4
)



plot_fancy_overlap_ppin <- function(gr, cancer, allele) {
    print(glue("plotting: {cancer} - {allele}"))
    p <- ggraph(gr, layout = "kk") +
        geom_edge_link(color = "grey75", width = 0.75, alpha = 0.7) +
        geom_node_point(
            aes(shape = interaction_source,
                color = node_color),
            size = 8
        ) +
        geom_node_text(aes(label = name), size = 3, family = "Arial") +
        scale_shape_manual(
            values = unlist(network_node_shapes),
            guide = guide_legend(title.position = "top",
                                 label.position = "top",
                                 order = 0)
        ) +
        scale_color_manual(
            values = network_node_colors,
            guide = guide_legend(title.position = "top",
                                 label.position = "top",
                                 order = 1)
        ) +
        theme_graph(base_family = "Arial", base_size = 8) +
        theme(
            legend.position = "bottom"
        ) +
        labs(
            title = glue("Interactions on the PPI for {allele} in {cancer}"),
            shape = "type of interaction",
            color = "details of interaction"
        )
    return(p)
}

fancy_overlap_ppin <- function(cancer, allele,
                               min_comp_size = 4) {
    set.seed(0)

    print(glue("beginning: {cancer} - {allele}"))

    synlet_df <- get_synthetic_lethal_data(cancer, allele) %>%
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
        filter(n() == which.max(abs(allele_diff))) %>%
        ungroup()
    genetic_df <- get_genetic_interaction_data(cancer, allele)

    df <- full_join(synlet_df, genetic_df, by = "hugo_symbol") %>%
        mutate(interaction_source = case_when(
            is.na(gene_cls) ~ "genetic",
            is.na(genetic_interaction) ~ "synthetic lethal",
            TRUE ~ "both"
        ))

    gr <- string_gr_prefilt %N>%
        filter(name %in% c(df$hugo_symbol, "KRAS")) %>%
        jhcutils::filter_component_size(min_size = min_comp_size)

    if (igraph::vcount(gr) == 0 | igraph::ecount(gr) == 0) return(NULL)

    save_graph_plot <- function(name, gr_plot, ...) {
        print(glue("saving: {cancer} - {allele}"))

        save_path <- plot_path(
            "40_10_overlap-synlet-comutation",
            paste0("overlap_ppi_", cancer, "_", allele, "_", name, ".svg")
        )
        ggsave_wrapper(gr_plot, save_path, width = 10, height = 8)
    }

    gr_components <- gr %N>%
        left_join(df, by = c("name" = "hugo_symbol")) %>%
        mutate(
            interaction_source = ifelse(
                is.na(interaction_source), "other", interaction_source
            ),
            synlet_direction = ifelse(
                allele_diff < 0, "synthetic lethal down", "synthetic lethal up"
            ),
            node_shape = network_node_shapes[interaction_source],
            node_color = case_when(
                interaction_source == "genetic" ~ genetic_interaction,
                interaction_source == "synthetic lethal" ~ synlet_direction,
                interaction_source == "other" ~ "other",
                interaction_source == "both" ~ "both"
            )
        ) %>%
        morph(to_components) %>%
        crystallize() %>%
        mutate(
            gr_plot = purrr::map(graph, plot_fancy_overlap_ppin,
                                 cancer = !!cancer, allele = !!allele)
        ) %>%
        pwalk(save_graph_plot)

}

depmap_gene_clusters_pairwise_df %>%
    select(cancer, group1, group2) %>%
    unique() %>%
    group_by(cancer) %>%
    summarise(allele = list(unique(c(unlist(group1), unlist(group2))))) %>%
    ungroup() %>%
    unnest(allele) %>%
    pwalk(fancy_overlap_ppin)
