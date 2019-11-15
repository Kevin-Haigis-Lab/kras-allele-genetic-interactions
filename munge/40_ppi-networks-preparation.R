################################
## Prepare STRING PPI network ##
################################

cache("string_gr", {

    info(logger, "Beginning preparation of the STRING PPI network.")

    string_data_path <- file.path("data", "STRING-network")

    edge_list_path <- file.path(string_data_path, "9606.protein.links.full.v11.0.txt")
    edge_list <- read_delim(edge_list_path, " ", n_max = Inf, col_types = cols(), progress = FALSE)
    string_gr <- as_tbl_graph(edge_list, directed = FALSE)

    name_tbl <- read_tsv(
        file.path(string_data_path, "human.name_2_string.tsv"),
        skip = 1,
        col_types = cols(),
        progress = FALSE,
        col_names = c("taxid", "hugo_symbol", "string_id")
    ) %>%
        select(-taxid)

    string_gr <- string_gr %N>%
        left_join(name_tbl, by = c("name" = "string_id")) %>%
        mutate(string_id = name,
               name = hugo_symbol) %>%
        select(-hugo_symbol)

    info(logger, glue("Number of nodes in STRING network: {igraph::vcount(string_gr)}"))
    info(logger, glue("Number of edges in STRING network: {igraph::ecount(string_gr)}"))
    info(logger, "Caching STRING network.")
    return(string_gr)
})
