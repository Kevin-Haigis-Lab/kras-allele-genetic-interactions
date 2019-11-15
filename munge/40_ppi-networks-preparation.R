##########################
## Prepare PPI networks ##
##########################

cache("string_gr",
{
    info(logger, "Beginning preparation of the STRING PPI network.")

    string_data_path <- file.path("data", "ppi-networks", "string")

    edge_list_path <- file.path(string_data_path,
                                "9606.protein.links.full.v11.0.txt")
    edge_list <- read_delim(edge_list_path, " ",
                            n_max = Inf,
                            col_types = cols(),
                            progress = FALSE)
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



cache("bioplex_gr",
{
    info(logger, "Beginning preparation of the BioPlex2 PPI network.")

    bioplex_data_path <- file.path("data", "ppi-networks", "bioplex2")

    edge_list_path <- file.path(bioplex_data_path, "bioplex2.txt")
    edge_list <- read_tsv(edge_list_path,
                            n_max = Inf,
                            col_types = cols(),
                            progress = FALSE) %>%
        janitor::clean_names() %>%
        mutate(from = symbol_a, to = symbol_b)

    bioplex_gr <- as_tbl_graph(edge_list, directed = FALSE)

    info(logger, glue("Number of nodes in BioPlex2 network: {igraph::vcount(bioplex_gr)}"))
    info(logger, glue("Number of edges in BioPlex2 network: {igraph::ecount(bioplex_gr)}"))
    info(logger, "Caching BioPlex2 network.")

    return(bioplex_gr)
})



# cache("intact_gr",
# {

# })



cache("hint_gr",
{
    info(logger, "Beginning preparation of the HINT PPI network.")

    hint_data_path <- file.path("data", "ppi-networks", "hint")

    edge_list_path_1 <- file.path(hint_data_path,
                                  "HomoSapiens_binary_hq.txt")
    edge_list_path_2 <- file.path(hint_data_path,
                                  "HomoSapiens_cocomp_hq.txt")
    edge_list <- bind_rows(
        {
            read_tsv(edge_list_path_1,
                     col_types = cols(),
                     progress = FALSE) %>% janitor::clean_names()
        },
        {
            read_tsv(edge_list_path_2,
                     col_types = cols(),
                     progress = FALSE) %>% janitor::clean_names()
        }
    ) %>%
        mutate(from = gene_a, to = gene_b)

    hint_gr <- as_tbl_graph(edge_list, directed = FALSE)

    info(logger, glue("Number of nodes in HINT network: {igraph::vcount(hint_gr)}"))
    info(logger, glue("Number of edges in HINT network: {igraph::ecount(hint_gr)}"))
    info(logger, "Caching HINT network.")

    return(hint_gr)
})
