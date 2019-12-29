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
prep_graph_for_plotting <- function(gr,
                                    geneset_tib,
                                    full_gs = NULL) {
    mod_gr <- gr %>%
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
    return(mod_gr)
}
prep_graph_for_plotting <- memoise::memoise(prep_graph_for_plotting)


# The main genes are colored by the interaction direction, and the
# edges' transparency is scaled by a combination of their centrality and if
# they are connected to one of the main genes.
plot_extracted_graph_of_geneset <- function(gr) {
    if (is.null(gr)) return(NULL)

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

    # Get names of the nodes on the shortest path.
    names(unlist(igraph::shortest_paths(gr, from, to, mode = "all")$vpath))
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
minimize_graph_to_shortest_paths <- memoise::memoise(minimize_graph_to_shortest_paths)


# Rank the nodes by their centrality.
# Make any nodes in `priority_nodes` infinite.
assign_centraility_rank <- function(gr, priority_nodes = NULL) {
    mod_gr <- gr %N>%
        mutate(
            node_ctrlty = centrality_pagerank(directed = FALSE),
            node_ctrlty = ifelse(name %in% !!priority_nodes, Inf, node_ctrlty),
            node_ctrlty = rank(node_ctrlty, ties.method = "random")
        )
    return(mod_gr)
}


# All node names in STRING.
STRING_NODE_NAMES <- igraph::V(string_gr)$name

# Extract the graph of the neighborhood of the genes driving a gene set's
# enrichment.
extract_local_neighborhood <- function(cancer, allele,
                                      datasource, geneset,
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
    # neighbors <- get_neighbors(enriched_geneset_tib$name, 2)

    gr <- minimize_graph_to_shortest_paths(
        string_gr, enriched_geneset_tib$name
    ) %>%
        get_giant_component() %>%
        assign_centraility_rank(priority_nodes = enriched_geneset_tib$name)

    max_ctrlty <- max(igraph::V(gr)$node_ctrlty)
    while (igraph::vcount(gr) > max_nodes) {

        # Remove the least central node.
        mod_gr <- gr %N>% filter(node_ctrlty != min(node_ctrlty))

        # If this splits the graph, fix and set centrality of removed node to
        # the maximum value.
        if (igraph::count_components(mod_gr) > 1) {
            gr <- gr %N>% mutate(node_ctrlty = ifelse(
                    node_ctrlty == min(node_ctrlty), !!max_ctrlty, node_ctrlty
                ))
        } else {
            gr <- mod_gr
        }

        # Break if all node centralities are the maximum value.
        if (all(igraph::V(gr)$node_ctrlty == max_ctrlty)) break
    }

    gr <- prep_graph_for_plotting(gr, enriched_geneset_tib,
                                  full_gs = all_geneset)

    return(gr)
}


# Clean the "glued" names for becoming a file name.
sanitize_save_names <- function(x) {
    str_replace_all(x, "\\/", "-") %>%
        str_replace_all(" ", "-") %>%
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
    save_name <- glue(glue_template)
    save_name <- sanitize_save_names(save_name)
    saveRDS(gr_plot, get_fig_proto_path(save_name, fig_num, supp = supp))
}



#### ---- Extract gene sets and their neighborhood from the STRING PPIN ---- ####

# A tibble of the gene sets to plot.
genesets_to_plot <- enrichr_tib %>%
    select(-gene_list) %>%
    unnest(cols = enrichr_res) %>%
    mutate(n_overlap = get_enrichr_overlap_int(overlap),
           term = standardize_enricher_terms(term)) %>%
    filter(adjusted_p_value < 0.05 & n_overlap >= 3) %>%
    select(cancer, allele, datasource, term) %>%
    unique() %>%
    dplyr::rename(geneset = term)

ProjectTemplate::cache("graphs_to_plot", depends = "enrichr_tib",
{
    NOT_NOVEL_GENES <- c("BRAF", "TP53", "APC", "NRAS", "PIK3CA",
                         "CTNNB1", "EGFR")
    graphs_to_plot <- pmap(genesets_to_plot,
                           extract_local_neighborhood,
                           max_nodes = 50,
                           min_interesting_genes = 3,
                           uninteresting_genes = NOT_NOVEL_GENES)
})

genesets_to_plot$gr <- graphs_to_plot
genesets_to_plot$gr_plot <- purrr::map(genesets_to_plot$gr,
                                       plot_extracted_graph_of_geneset)

save_tmplt <- "{cancer}_{allele}_{datasource}_{geneset}"
genesets_to_plot %>%
    filter(!is.null(gr_plot)) %>%
    pwalk(save_neighborhood_plot, glue_template = save_tmplt)


#### ---- Select specific examples for saving as protos for figures. ---- ####

protos_save_tib <- tibble::tribble(
    ~cancer, ~allele, ~datasource, ~geneset, ~fig_num, ~supp,
    "LUAD", "G12C", "Transcription_Factor_PPIs", "MYC", 11, TRUE,
    "LUAD", "G12D", "KEGG_2019_Human", "Focal adhesion", 11, TRUE
)


# Make any adjustments to the graph specifically for the plots used
# for figures.
adjust_gr_for_figure_proto <- function(gr) {
    mod_gr <- gr %N>%
        mutate(
            node_size = node_size / 2,
            label_size = label_size / 2
        )
    return(mod_gr)
}

genesets_to_plot %>%
    filter(!is.null(gr_plot)) %>%
    right_join(protos_save_tib,
               by = c("cancer", "allele", "datasource", "geneset")) %>%
    mutate(
        gr = purrr::map(gr, adjust_gr_for_figure_proto),
        gr_plot = purrr::map(gr, plot_extracted_graph_of_geneset)
    ) %>%
    pwalk(save_neighborhood_proto, glue_template = save_tmplt)

