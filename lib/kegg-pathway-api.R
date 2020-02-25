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
        mutate(pretty_name = ifelse(!is.na(hugo_symbol),
                                  hugo_symbol,
                                  display_name)) %>%
        select(pretty_name, everything())
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


node_uniqueid_mapping <- function(nodes_df) {
    node_idx <- nodes_df %>%
        select(pretty_name) %>%
        unique() %>%
        arrange(pretty_name) %>%
        mutate(unique_id = seq(1, n()))

    nodes_df %>%
        select(pretty_name, entry_id) %>%
        unique() %>%
        left_join(node_idx, by = "pretty_name") %>%
        select(unique_id, entry_id)
}


replace_entryid_with_uniqueid <- function(edges_df, unique_map) {
    replace_ids <- function(eids, uid_map) {
        tibble(entry_id = eids) %>%
            left_join(uid_map, by = "entry_id") %>%
            pull(unique_id)
    }

    new_edge_df <- edges_df %>%
        mutate(entry_id_from = from,
               entry_id_to = to,
               from = replace_ids(from, unique_map),
               to = replace_ids(to, unique_map))
    return(new_edge_df)
}



# Parse the KGML downloaded from KEGG Pathways into a tidygraph.
parse_kegg_kgml <- function(path) {
    kgr <- KEGGgraph::parseKGML(path)
    nodes <- parse_kegg_nodes(KEGGgraph::nodes(kgr))
    edges <- parse_kegg_edges(KEGGgraph::edges(kgr))

    entryid_uniqueid_map <- node_uniqueid_mapping(nodes)
    edges <- replace_entryid_with_uniqueid(edges, entryid_uniqueid_map)
    nodes <- left_join(nodes, entryid_uniqueid_map, by = "entry_id") %>%
        mutate(name = unique_id) %>%
        group_by(name, pretty_name, node_type) %>%
        summarise(entry_ids = list(entry_id)) %>%
        ungroup()
        
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


# An opinionated way to process and output a KEGG pathway.
parse_and_annotate_kegg_kgml <- function(path) {
    gr <- parse_kegg_kgml(path) %N>%
        mutate(name = pretty_name) %>%
        select(-pretty_name) %>%
        filter(!str_detect(name, "TITLE") &
               !(node_type %in% c("map", "compound"))) %>%
        annotate_kegg_edges() %>%
        filter(interaction_subtype_name != "compound")
    return(gr)
}
