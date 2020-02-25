# Take a KEGG pathway and convert to a tidygraph object.

# Convert KEGG ID ("hsa:####") to Hugo Symbol.
keggid_to_hugosymbol <- function(kegg_id) {
    mget(KEGGgraph::translateKEGGID2GeneID(kegg_id),
         org.Hs.egSYMBOL,
         ifnotfound = NA_character_)
}


# Parse a node from a `KEGGpathway` object from 'KEGGgraph'.
# Returns a tibble with the information.
parse_kegg_node <- function(node) {
    tibble(
        name = node@entryID,
        kegg_name = node@name,
        hugo_name = unlist(keggid_to_hugosymbol(node@name)),
        display_name = node@graphics@name,
        node_type = node@type,
        reaction = node@reaction,
        node_map = node@map,
        graphics_x = node@graphics@x,
        graphics_y = node@graphics@y,
        graphics_shape = node@graphics@type,
        graphics_width = node@graphics@width,
        graphics_height = node@graphics@height
    )
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
    purrr::map(kegg_nodes, parse_kegg_node) %>%
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


# Parse the KGML downloaded from KEGG Pathways.
parse_kegg_kgml <- function(path) {
    kgr <- KEGGgraph::parseKGML(path)
    nodes <- parse_kegg_nodes(KEGGgraph::nodes(kgr))
    edges <- parse_kegg_edges(KEGGgraph::edges(kgr))

    tbl_graph(
        nodes = nodes,
        edges = edges
    )
}


# Annotate the edges with the subtype interaction.
annotate_kegg_edges <- function(gr) {
    gr %E>% mutate(
        interaction_subtype_name = map_chr(interaction_subtype, ~ .x@name),
        interaction_subtype_value = map_chr(interaction_subtype, ~ .x@value)
    )
}


# test_fpath <- file.path(
#     "data", "kegg-pathways", "hsa04010_MAPK-signaling-pathway.xml"
# )

# parse_kegg_kgml(test_fpath) %N>%
#     mutate(name = ifelse(is.na(hugo_name), display_name, hugo_name)) %>%
#     annotate_kegg_edges()
