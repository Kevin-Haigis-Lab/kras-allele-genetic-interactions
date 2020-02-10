# Make a data structure for the overlap of the synthetic lethal and genetic
# interaction analyses.

GRAPHS_DIR <- "40_10_overlap-synlet-comutation"
reset_graph_directory(GRAPHS_DIR)

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



#### ---- Plotting ---- ####

# Dictionary of node shapes to use by genetic interaction.
network_node_shapes <- list(
    "comutation" = 15,
    "dependency" = 17,
    "both" = 18,
    "other" = 19
)

# Dictionary of node colors based on direction of genetic interaction.
network_node_colors <- list(
    "increased\ncomut." = comut_mutex_pal[["comutation"]],
    "reduced\ncomut." = comut_mutex_pal[["exclusivity"]],
    "reduced dep." = synthetic_lethal_pal[["up"]],
    "increased dep." = synthetic_lethal_pal[["down"]],
    "both" = "#ED6165",
    "other" = "#F44C4C"
) %>%
    lighten(factor = 1.4)

# Make the purple even brighter.
network_node_colors["increased dep."] <- lighten(
    network_node_colors[["increased dep."]],
    factor = 1.4
)


# The simple graph of the combined graph.
# The graph is the combination of STRING, HINT, and BioPlex2.
simple_combined_ppi_gr <- convert(combined_ppi_gr, to_simple) %E>%
    mutate(num_source = purrr::map_dbl(.orig_data,
                                       ~ n_distinct(.x$source))) %>%
    select(-.tidygraph_edge_index, -.orig_data) %N>%
    select(-.tidygraph_node_index)


# Plot a graph with annotations for comutation and genetic dependencies.
plot_fancy_overlap_ppin <- function(gr, cancer, allele) {
    print(glue("plotting: {cancer} - {allele}"))
    p <- ggraph(gr, layout = "kk") +
        geom_edge_link(
            aes(color = num_source,
                width = num_source),
            alpha = 0.7) +
        scale_edge_width_continuous(
            range = c(0.75, 1.5),
            guide = guide_legend(title.position = "top",
                                 label.position = "top",
                                 order = 3)) +
        scale_edge_color_continuous(
            low = "gray80", high = "gray50",
            guide = FALSE) +
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
                                 order = 1)
        ) +
        scale_color_manual(
            values = network_node_colors,
            guide = guide_legend(title.position = "top",
                                 label.position = "top",
                                 order = 2)
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


# Assign the node shape and color based on a combination of the comutation
# and the dependency analysis.
assign_node_color_and_shape <- function(gr) {
    new_gr <- gr %N>%
        mutate(
            genetic_interaction = ifelse(
                is.na(genetic_interaction), NA_character_, genetic_interaction
            ),
            interaction_source = ifelse(
                is.na(interaction_source), "other", interaction_source
            ),
            synlet_direction = ifelse(
                allele_diff < 0, "increased dep.", "reduced dep."
            ),
            node_shape = network_node_shapes[interaction_source],
            node_color = case_when(
                interaction_source == "comutation" ~ genetic_interaction,
                interaction_source == "dependency" ~ synlet_direction,
                interaction_source == "other" ~ "other",
                interaction_source == "both" ~ "both"
            )
        )
    return(new_gr)
}


fancy_overlap_ppin <- function(cancer, allele,
                               min_comp_size = 4,
                               ignore_genes = c()) {
    set.seed(0)

    print(glue("beginning: {cancer} - {allele}"))

    res <- get_overlapped_gr(cancer, allele, min_comp_size, ignore_genes)
    gr <- res[["graph"]]
    df <- res[["data"]]

    if (igraph::vcount(gr) == 0 | igraph::ecount(gr) == 0) return(NULL)

    save_graph_plot <- function(name, gr_plot, ...) {
        print(glue("saving: {cancer} - {allele}"))

        save_path <- plot_path(
            GRAPHS_DIR,
            paste0("overlap_ppi_", cancer, "_", allele, "_", name, ".svg")
        )
        ggsave_wrapper(gr_plot, save_path, width = 10, height = 8)
    }

    gr_components <- gr %N>%
        left_join(df, by = c("name" = "hugo_symbol")) %>%
        assign_node_color_and_shape() %>%
        morph(to_components) %>%
        crystallize() %>%
        mutate(
            gr_plot = purrr::map(graph, plot_fancy_overlap_ppin,
                                 cancer = !!cancer, allele = !!allele)
        ) %>%
        pwalk(save_graph_plot)

}


clustered_fancy_overlap_ppin <- function(cancer, allele,
                               min_comp_size = 4,
                               ignore_genes = c()) {
    set.seed(0)
    print(glue("beginning (clustered): {cancer} - {allele}"))

    res <- get_overlapped_gr(cancer, allele, min_comp_size, ignore_genes)
    gr <- res[["graph"]]
    df <- res[["data"]]

    if (igraph::vcount(gr) == 0 | igraph::ecount(gr) == 0) return(NULL)

    save_graph_plot <- function(name, gr_plot, ...) {
        print(glue("saving (clustered): {cancer} - {allele}"))

        save_path <- plot_path(
            GRAPHS_DIR,
            paste0("clustered_ppi_", cancer, "_", allele, "_", name, ".svg")
        )
        ggsave_wrapper(gr_plot, save_path, width = 10, height = 8)
    }

    gr_components <- gr %N>%
        left_join(df, by = c("name" = "hugo_symbol")) %>%
        assign_node_color_and_shape() %N>%
        morph(to_components) %>%
        mutate(cls = group_spinglass()) %>%
        unmorph() %E>%
        filter(.N()$cls[from] == .N()$cls[to]) %N>%
        morph(to_components) %>%
        crystallize() %>%
        mutate(
            gr_plot = purrr::map(graph, plot_fancy_overlap_ppin,
                                 cancer = !!cancer, allele = !!allele)
        ) %>%
        pwalk(save_graph_plot)

}


genes_to_ignore <- c("TTN")

depmap_gene_clusters_pairwise_df %>%
    select(cancer, group1, group2) %>%
    unique() %>%
    group_by(cancer) %>%
    summarise(allele = list(unique(c(unlist(group1), unlist(group2))))) %>%
    ungroup() %>%
    unnest(allele) %>%
    pwalk(fancy_overlap_ppin, ignore_genes = genes_to_ignore) %>%
    pwalk(clustered_fancy_overlap_ppin, ignore_genes = genes_to_ignore)




#### ---- Overlap of comutation and genetic dependency analysis ---- ####

get_shared_comutation_dependency <- function(cancer, allele) {
    get_overlapped_df(cancer, allele) %>%
        select(-c(group1, group2, g1_avg, g2_avg,
                  g1_other_avg, g2_other_avg)) %>%
        filter(interaction_source == "both") %>%
        mutate(cancer = !!cancer) %>%
        select(cancer, allele, everything())
}


depmap_gene_clusters_pairwise_df %>%
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


#### ---- Comparing of COAD allele overlap networks ---- ####




annotate_merged_graph <- function(gr, df) {
    for (i in 1:nrow(df)) {
        allele_i <- df$allele[[i]]
        nodes <- igraph::V(df$ppi[[i]])$name
        gr <- gr %N>%
            mutate(allele = case_when(
                name %in% !!nodes & allele == "" ~ !!allele_i,
                name %in% !!nodes ~ paste(allele, !!allele_i, sep = ", "),
                TRUE ~ allele
            ))
    }
    return(gr)
}


make_gray_pal <- function(max_ct, start = 0.8, end = 0.5) {
    pal <- grDevices::gray.colors(max_ct - 1, start = start, end = end)
    names(pal) <- paste(seq(2, max_ct), "alleles", sep = " ")
    return(pal)
}


factor_multiple_node_colors <- function(nc) {
    alleles <- unique(nc[!str_detect(nc, "allele")])
    grays <- unique(nc[str_detect(nc, "allele")])
    levels <- c(sort(alleles), sort(grays))
    return(factor(nc, levels = levels))
}


get_allele_count_label <- function(allele) {
    paste(str_count(allele, ",") + 1, "alleles", sep = " ")
}


plot_overlap_comparison_graph <- function(gr, special_labels = NULL) {
    max_ct <- max(str_count(igraph::V(gr)$allele, ",")) + 1
    pal <- c(short_allele_pal,
             make_gray_pal(max_ct),
             "special" = "indianred1")
    p <- gr %N>%
        mutate(
            node_color = case_when(
                name %in% !!special_labels ~ "special",
                str_detect(allele, ",") ~ get_allele_count_label(allele),
                TRUE ~ allele
            ),
            node_color = factor_multiple_node_colors(node_color)
        ) %>%
        ggraph(layout = "kk") +
        geom_edge_link(width = 0.3, color = "gray30", alpha = 0.4) +
        geom_node_point(aes(color = node_color), size = 3) +
        geom_node_text(aes(label = name), size = 2, family = "Arial", repel = TRUE) +
        scale_color_manual(values = pal) +
        theme_graph()
    return(p)
}


make_overlap_comparison_graph <- function(df) {
    merged_gr <- recursive_graph_join(df$ppi) %>%
        convert(to_simple, .clean = TRUE) %>%
        convert(to_undirected, .clean = TRUE) %N>%
        select(name) %E>%
        select(from, to) %N>%
        mutate(allele = "") %>%
        annotate_merged_graph(df)
    return(merged_gr)
}


special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53")

coad_overlap_comparison_plot <- tibble(cancer = c("COAD", "COAD", "COAD"),
       allele = c("G12D", "G12V", "G13D")) %>%
    mutate(
        ppi = purrr::map2(cancer, allele,
                          get_overlapped_gr,
                          min_comp_size = 4, ignore_genes = genes_to_ignore),
        data = purrr::map(ppi, ~ .x$data),
        ppi = purrr::map(ppi, ~.x$graph)
    ) %>%
    make_overlap_comparison_graph() %>%
    plot_overlap_comparison_graph(special_labels = special_nodes)
ggsave_wrapper(
    coad_overlap_comparison_plot,
    plot_path(GRAPHS_DIR, "coad_overlap_comparison_plot.svg"),
    "large"
)
