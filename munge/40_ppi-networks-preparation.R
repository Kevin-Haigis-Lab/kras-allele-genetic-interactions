##########################
## Prepare PPI networks ##
##########################

ProjectTemplate::cache("string_gr",
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

    string_gr <- string_gr  %E>%
        filter(combined_score > 800) %N>%
        left_join(name_tbl, by = c("name" = "string_id")) %>%
        mutate(string_id = name,
               name = hugo_symbol) %>%
        select(-hugo_symbol)

    info(logger, glue("Number of nodes in STRING network: {igraph::vcount(string_gr)}"))
    info(logger, glue("Number of edges in STRING network: {igraph::ecount(string_gr)}"))
    info(logger, "Caching STRING network.")
    return(string_gr)
})



ProjectTemplate::cache("bioplex_gr",
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



# ProjectTemplate::cache("intact_gr",
# {

# })



ProjectTemplate::cache("hint_gr",
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



# Some summary stats of for the graphs.
ProjectTemplate::cache("ppi_graph_summary_stats",
                       depends = c("string_gr", "bioplex_gr", "hint_gr"),
{
    # Calculates, prints, and returns (as a tibble) summary statistics
    #   for a graph (igraph or tidygraph objects.
    print_graph_summary_stats <- function(gr, name) {

        num_edges <- igraph::vcount(gr)
        num_vertices <- igraph::ecount(gr)
        num_mult <- sum(igraph::count_multiple(gr) > 1)
        gr_diameter <- igraph::diameter(gr)
        gr_clqs <- igraph::count_max_cliques(gr)
        gr_dist <- igraph::mean_distance(gr)
        gr_trans <- igraph::transitivity(gr, type = 'global')

        cat("Summary stats for", name, "\n")
        cat("  number of nodes:", num_edges, "\n")
        cat("  number of edges:", num_vertices, "\n")
        cat("  number of multiple or loop edges:", num_mult, "\n")
        cat("  diameter:", gr_diameter, "\n")
        cat("  num. maximum cliques:", gr_clqs, "\n")
        cat("  mean distance:", gr_dist, "\n")
        cat("  global transitivity:", gr_trans, "\n")
        cat("---------\n")

        tibble(
            graph_name = name,
            metric_name = c("num_edges", "num_vertices", "num_multiple",
                            "diameter", "num_max_cliques", "mean_distance",
                            "transitivity"),
            value = c(num_edges, num_vertices, num_mult, gr_diameter,
                      gr_clqs, gr_dist, gr_trans)
        )
    }

    ppi_graph_summary_stats <- purrr::map2(
        list(string_gr, bioplex_gr, hint_gr),
        list("STRING", "BioPlex2", "HINT"),
        print_graph_summary_stats
    ) %>%
        bind_rows()

    return(ppi_graph_summary_stats)
})

ppi_graph_summary_stats %>%
    pivot_wider(names_from = graph_name, values_from = value)



ProjectTemplate::cache("combined_ppi_gr",
                       depends = c("string_gr", "bioplex_gr", "hint_gr"),
{
    combined_ppi_gr <- graph_join(
        {
            string_gr %E>%
                select(from, to, combined_score) %>%
                mutate(source = "STRING")
        },
        {
            bioplex_gr %E>%
                select(from, to, p_interaction) %>%
                mutate(source = "BioPlex2")
        },
        by = "name"
    ) %N>%
        graph_join(
            { hint_gr %E>% select(from, to) %>% mutate(source = "HINT") },
            by = "name"
        )

    return(combined_ppi_gr)
})
