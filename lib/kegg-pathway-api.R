# Take a KEGG pathway and convert to a tidygraph object.

# Convert KEGG ID ("hsa:####") to Hugo Symbol.
keggid_to_hugosymbol <- function(kegg_id) {
    mget(KEGGgraph::translateKEGGID2GeneID(kegg_id),
         org.Hs.egSYMBOL,
         ifnotfound = NA_character_)
}


# Add a column `name` to the nodes data frame.
add_pretty_node_name <- function(nodes_df) {
    nodes_df %>%
        mutate(name = ifelse(!is.na(hugo_symbol),
                                  hugo_symbol,
                                  display_name)) %>%
        select(name, everything())
}


# Parse a node from a `KEGGpathway` object from 'KEGGgraph'.
# Returns a tibble with the information.
parse_kegg_node <- function(node) {
    tibble(
        hugo_symbol = unlist(keggid_to_hugosymbol(node@name))[[1]],
        entry_id = node@entryID,
        kegg_name = node@name,
        display_name = node@graphics@name,
        node_type = node@type,
        reaction = node@reaction,
        node_map = node@map,
        graphics_x = node@graphics@x,
        graphics_y = node@graphics@y,
        graphics_shape = node@graphics@type,
        graphics_width = node@graphics@width,
        graphics_height = node@graphics@height
    ) %>%
        add_pretty_node_name()
}


parse_kegg_node_simple <- function(node) {
    tibble(
        hugo_symbol = unlist(keggid_to_hugosymbol(node@name)),
        entry_id = node@entryID,
        kegg_name = node@name,
        display_name = node@graphics@name,
        node_type = node@type
    ) %>%
        add_pretty_node_name()
}


# Parse an edge from a `KEGGpathway` object from 'KEGGgraph'.
# Returns a tibble with the information.
parse_kegg_edge <- function(edge) {
    tibble(
        from = edge@entry1ID,
        to = edge@entry2ID,
        interaction_type = edge@type,
        interaction_subtype = edge@subtype
    )
}


# Parse the nodes from a `KEGGpathway` object from 'KEGGgraph'.
# Returns a tibble with the information.
parse_kegg_nodes <- function(kegg_nodes) {
    purrr::map(kegg_nodes, parse_kegg_node_simple) %>%
        bind_rows() %>%
        unique()
}


# Parse the edges from a `KEGGpathway` object from 'KEGGgraph'.
# Returns a tibble with the information.
parse_kegg_edges <- function(kegg_edges) {
    purrr::map(kegg_edges, parse_kegg_edge) %>%
        bind_rows() %>%
        unique()
}


get_relevent_edges <- function(name, entry_id, edge_df, col) {
    col <- rlang::enquo(col)
    edge_df %>%
        filter(!!col == !!entry_id) %>%
        mutate(!!col := !!name)
}


# Parse the KGML downloaded from KEGG Pathways into a tidygraph.
parse_kegg_kgml <- function(path) {
    kgr <- KEGGgraph::parseKGML(path)
    nodes <- parse_kegg_nodes(KEGGgraph::nodes(kgr))
    edges <- parse_kegg_edges(KEGGgraph::edges(kgr))

    mod_edges <- nodes %>%
        select(name, entry_id) %>%
        unique() %>%
        pmap(get_relevent_edges, edges, from) %>%
        bind_rows() %>%
        unique()
    mod_edges <- nodes %>%
        select(name, entry_id) %>%
        unique() %>%
        pmap(get_relevent_edges, mod_edges, to) %>%
        bind_rows() %>%
        unique()

    nodes <- nodes %>%
        select(name, node_type) %>%
        unique()

    tbl_graph(
        nodes = nodes,
        edges = mod_edges
    )
}


# Annotate the edges with the subtype interaction.
annotate_kegg_edges <- function(gr) {
    gr %E>% mutate(
        interaction_subtype_name = map_chr(interaction_subtype, ~ .x@name),
        interaction_subtype_value = map_chr(interaction_subtype, ~ .x@value)
    )
}


# An opinionated way to process and output a KEGG pathway.
parse_and_annotate_kegg_kgml <- function(path) {
    gr <- parse_kegg_kgml(path) %N>%
        filter(!str_detect(name, "TITLE") &
               !(node_type %in% c("map", "compound"))) %>%
        annotate_kegg_edges() %>%
        filter(interaction_subtype_name != "compound")
    return(gr)
}
