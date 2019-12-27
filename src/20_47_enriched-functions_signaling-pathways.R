# Explicit and specific signaling pathways from the functional enrichment.

GRAPHS_DIR <- "20_47_enriched-functions_signaling-pathways"
reset_graph_directory(GRAPHS_DIR)

library(ggraph)

#### ---- General Setup ---- ####

enrichr_tib_mod <- enrichr_tib %>%
    select(-gene_list) %>%
    unnest(cols = enrichr_res) %>%
    mutate(n_overlap = get_enrichr_overlap_int(overlap)) %>%
    mutate(overlap_genes = str_split(genes, ";")) %>%
    select(-genes) %>%
    mutate(term_mod = standardize_enricher_terms(term))


#### ---- Subroutines ---- ####

# Get the genes underlying the enrichment of a gene set.
get_overlap_genes <- function(cancers, alleles = NULL, geneset_names) {
    enrichr_tib_mod %>%
        filter(cancer %in% !!cancers & term_mod %in% !!geneset_names) %>%
        filter(is.null(!!alleles) | allele %in% alleles) %>%
        pull(overlap_genes) %>%
        unlist() %>%
        unique()
}


# Return the type of interaction between an allele and other gene.
# Returns either "increased", "reduced", or NULL.
get_interaciton_type <- function(cancer, allele, gene, stop_if_empty = TRUE) {
    interaction_type <- genetic_interaction_df %>%
        filter(
            cancer == !!cancer & allele == !!allele & hugo_symbol == !!gene
        ) %>%
        pull(genetic_interaction) %>%
        switch_comut_terms()
    if (stop_if_empty & length(interaction_type) == 0) {
        msg <- "No interaction type found for {allele} & {gene} in {cancer}."
        stop(glue(msg))
    }
    return(interaction_type)
}


# Vectorize `get_interaciton_type()` over the `genes`.
get_interaciton_types <- function(genes, cancer, allele, stop_if_empty = TRUE) {
    purrr::map_chr(genes, get_interaciton_type,
                   cancer = cancer,
                   allele = allele,
                   stop_if_empty = stop_if_empty)
}


# Get the neighbors for `genes` in STRING.
get_neighbors <- function(genes, min_connections = 2) {
    string_gr %E>%
        filter(
            .N()$name[from] %in% !!genes |
            .N()$name[to] %in% !!genes
        ) %N>%
        filter(centrality_degree(mode = "all") >= !!min_connections) %>%
        as_tibble() %>%
        u_pull(name)
}
get_neighbors <- memoise::memoise(get_neighbors)


# Get the genes for a term from a data source.
geneset_genes <- function(datasource, term) {
    genes <- enrichr_genesets %>%
        filter(datasource == !!datasource & std_term == !!term)
    return(unname(unlist(genes)))
}
geneset_genes <- memoise::memoise(geneset_genes)


# Extract the subnetwork from STRING and add attributes for plotting.
extract_graph_for_plotting <- function(geneset_tib,
                                       neighbors,
                                       full_gs = NULL) {
    gr <- string_gr %N>%
        filter(name %in% !!neighbors | name %in% !!geneset_tib$name) %>%
        left_join(geneset_tib, by = "name") %>%
        mutate(
            interaction = case_when(
                !is.na(interaction) ~ interaction,
                name %in% full_gs ~ "in_geneset",
                TRUE ~ "none"
            ),
            node_size = ifelse(name %in% geneset_tib$name, 4, 2),
            label_size = ifelse(name %in% geneset_tib$name, 4, 3),
            label_face = ifelse(interaction == "none", "plain", "bold")
        ) %E>%
        mutate(
            edge_centrality = centrality_edge_betweenness(directed = FALSE),
            edge_centrality = scales::rescale(edge_centrality, to = c(0, 1)),
            edge_is_link1 = .N()$name[from] %in% geneset_tib$name,
            edge_is_link2 = .N()$name[to] %in% geneset_tib$name,
            edge_is_link = edge_is_link1 + edge_is_link2,
            edge_alpha = (edge_is_link2 * 2) + edge_centrality,
            edge_alpha = scales::rescale(edge_alpha, to = c(0.1, 0.7))
        )
    return(gr)
}
extract_graph_for_plotting <- memoise::memoise(extract_graph_for_plotting)


# Plot the graph of the extracted gene set and neighbors.
plot_extracted_graph_of_geneset <- function(gr) {
    p <- ggraph(gr, layout = "stress") +
        geom_edge_link(
            aes(alpha = edge_alpha),
            color = "grey50",
        ) +
        scale_edge_alpha_identity() +
        geom_node_point(
            aes(color = interaction,
                size = node_size)
        ) +
        geom_node_text(
            aes(label = name,
                fontface = label_face,
                size = label_size),
            repel = TRUE,
            family = "Arial") +
        scale_color_manual(
            values = c(comut_updown_pal,
                       "none" = "grey70",
                       "in_geneset" = "grey40")
        ) +
        scale_size_identity() +
        theme_graph(base_size = 7, base_family = "Arial")
    return(p)
}


# Get all pairs of the nodes in `gs`.
get_all_node_pairs <- function(gs) {
    combn(gs, 2) %>%
        t() %>%
        as.data.frame(strings_as_factors = FALSE) %>%
        as_tibble() %>%
        set_colnames(c("from", "to"))
}


# Returns a boolean as to whether a path exists between `from` and `to`.
# `from` and `to` must be node indices.
path_does_exist <- function(gr, from, to) {
    !any(is.infinite(igraph::distances(gr, from, to)))
}


# Get the nodes along the shortest path(s) between `from` and `to`.
# Path the `from` and `to` as node names.
get_nodes_in_shortest_path <- function(from, to, gr) {
    from_idx <- get_node_index(gr, name == from)
    to_idx <- get_node_index(gr, name == to)

    # Skip if there are no shortest paths.
    if (!path_does_exist(gr, from_idx, to_idx)) return(c(from, to))
    gr %>%
        convert(to_shortest_path,
                from = from_idx, to = to_idx,
                mode = "all") %N>%
        as_tibble() %>%
        pull(name)
}


# Shrink the graph to only include the nodes that lie on the geodesic path
# of the main nodes.
minimize_graph_to_shortest_paths <- function(gr, main_nodes) {
    node_pairs <- get_all_node_pairs(main_nodes)
    genes_on_paths <- pmap(node_pairs, get_nodes_in_shortest_path, gr = gr) %>%
        unlist() %>%
        unique()
    mod_gr <- gr %N>%
        filter(name %in% genes_on_paths)
    return(mod_gr)
}


# All node names in STRING.
STRING_NODE_NAMES <- as_tibble(string_gr, active = 'nodes')$name

# Plot the graph of the neighborhood of the genes driving a gene set's
# enrichment. The main genes are colored by the interaction direction, and the
# edges' transparency is scaled by a combination of their centrality and if
# they are connected to one of the main genes.
plot_labeled_neighborhood <- function(cancer, allele,
                                      datasource, geneset,
                                      max_neighbors = Inf,
                                      min_interesting_genes = 0,
                                      uninteresting_genes = NULL,
                                      max_nodes = 100,
                                      ...) {
    set.seed(0)  # For reproducibility

    # The genes and the comutation interaction that are enriched in the geneset.
    enriched_geneset_tib <- tibble(
            name = get_overlap_genes(cancer, allele, geneset)
        ) %>%
        mutate(interaction = get_interaciton_types(name, cancer, allele)) %>%
        filter(name %in% !!STRING_NODE_NAMES)

    n_interesting_genes <- sum(!enriched_geneset_tib$name %in% uninteresting_genes)
    if (n_interesting_genes < min_interesting_genes) {
        return(NULL)
    }

    # All genes in gene set.
    all_geneset <- geneset_genes(datasource, geneset)

    # Neighbors of the genes.
    neighbors <- get_neighbors(enriched_geneset_tib$name, 2)
    if (length(neighbors) > max_neighbors) {
        neighbors <- get_neighbors(enriched_geneset_tib$name, 3)
    }

    # Get the graph with all neighbors of the enriched genes and all edges.
    gr <- extract_graph_for_plotting(enriched_geneset_tib,
                                     neighbors,
                                     all_geneset)

    if (igraph::vcount(gr) > max_nodes) {
        gr <- minimize_graph_to_shortest_paths(gr, enriched_geneset_tib$name)
    }


    if (igraph::vcount(gr) > max_nodes) {
        gr <- gr %N>%
            filter(name %in% c(enriched_geneset_tib$name, all_geneset))
    }

    # Plot the graph.
    gr_plot <- plot_extracted_graph_of_geneset(gr)
    return(gr_plot)
}


# Clean the "glued" names for becoming a file name.
sanitize_save_names <- function(x) {
    str_replace_all(x, "\\/", "-") %>%
        str_remove_all("\\.")
}


# The the plot as an image.
save_neighborhood_plot <- function(gr_plot,
                                   cancer,
                                   allele,
                                   datasource,
                                   geneset,
                                   glue_template,
                                   ...) {
    if (is.null(gr_plot)) return(NULL)

    datasource <- str_replace_all(datasource, "_", "-")
    geneset <- str_replace_all(geneset, ' ', '-')
    save_name <- glue(glue_template)
    save_name <- sanitize_save_names(save_name)
    save_name <- paste0(save_name, ".svg")
    ggsave_wrapper(gr_plot, plot_path(GRAPHS_DIR, save_name), "large")
}


# Save the plot as a proto for a figure.
save_neighborhood_proto <- function(gr_plot,
                                   cancer,
                                   allele,
                                   datasource,
                                   geneset,
                                   fig_num,
                                   supp,
                                   glue_template,
                                   ...) {
    if (is.null(gr_plot)) return(NULL)

    datasource <- str_replace_all(datasource, "_", "-")
    geneset <- str_replace_all(geneset, ' ', '-')
    save_name <- glue(glue_template)
    save_name <- sanitize_save_names(save_name)
    saveRDS(gr_plot, get_fig_proto_path(save_name, fig_num, supp = supp))
}



#### ---- Extract gene sets and their neighborhood from the STRING PPIN ---- ####

NOT_NOVEL_GENES <- c("BRAF", "TP53", "APC", "NRAS", "PIK3CA", "CTNNB1", "EGFR")

genesets_to_plot <- enrichr_tib %>%
    select(-gene_list) %>%
    unnest(cols = enrichr_res) %>%
    mutate(n_overlap = get_enrichr_overlap_int(overlap),
           term = standardize_enricher_terms(term)) %>%
    filter(adjusted_p_value < 0.05 & n_overlap >= 3) %>%
    select(cancer, allele, datasource, term) %>%
    unique() %>%
    dplyr::rename(geneset = term) %>%
    filter(cancer == "LUAD")

genesets_to_plot$gr_plot <- pmap(genesets_to_plot, plot_labeled_neighborhood,
                                 min_interesting_genes = 3,
                                 uninteresting_genes = NOT_NOVEL_GENES)

save_tmplt <- "{cancer}_{allele}_{datasource}_{geneset}"
pwalk(genesets_to_plot, save_neighborhood_plot, glue_template = save_tmplt)


# TODO: Select specific examples for saving as protos for figures.
# pwalk(genesets_to_plot, save_neighborhood_proto, glue_template = save_tmplt)
