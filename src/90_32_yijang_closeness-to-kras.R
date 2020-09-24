# Yi-jang asked how close a set of genes were to KRas on the PPIN

GRAPHS_DIR <- "90_32_yijang_closeness-to-kras"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


set.seed(0)

YJ_GENES <- c("HSPA13", "TP53BP2", "ZKSCAN8")

if (all(YJ_GENES %in% igraph::V(combined_ppi_gr)$name)) {
  message("All of Y-Jang's gene's are in the combined PPIN.")
} else {
  stop("Not all of Y-Jang's gene's are in the combined PPIN.")
}


match_length <- function(x, y) {
  if (length(x) == 1 & length(y) > 1) {
    x <- rep(x, length(y))
  }
  return(x)
}


shortest_distance <- function(gr, n1, n2) {
  igraph::distances(gr, n1, n2)[1, 1]
}


shortest_distances <- function(gr, n1, n2) {
  n1 <- match_length(n1, n2)
  map2_dbl(n1, n2, ~ shortest_distance(gr, .x, .y))
}


get_giant_component_memo <- memoise::memoise(get_giant_component)

random_shortest_distances <- function(gr, n1, n = 100) {
  f <- function(x, gr) {
    random_nodes <- sample(igraph::V(get_giant_component_memo(gr))$name, n)
    y <- shortest_distances(gr, x, random_nodes)
    return(y[is.finite(y)])
  }
  map(n1, f, gr = gr)
}


number_shortest_path <- function(gr, n1, n2) {
  length(igraph::all_shortest_paths(combined_ppi_gr, n1, n2)$res)
}


number_shortest_paths <- function(gr, n1, n2) {
  n1 <- match_length(n1, n2)
  map2_dbl(n1, n2, ~ number_shortest_path(gr, .x, .y))
}


random_number_shortest_paths <- function(gr, n1, n = 100) {
  f <- function(x, gr) {
    random_nodes <- sample(igraph::V(get_giant_component_memo(gr))$name, n)
    y <- number_shortest_paths(gr, x, random_nodes)
    return(y[is.finite(y)])
  }
  map(n1, f, gr = gr)
}


extract_shortest_path <- function(gr, n1, n2) {
  geodesic_nodes <- igraph::all_shortest_paths(combined_ppi_gr, n1, n2)$res %>%
    unlist() %>%
    names() %>%
    unique()
  gr %N>%
    filter(name %in% !!geodesic_nodes)
}


extract_shortest_paths <- function(gr, n1, n2) {
  n1 <- match_length(n1, n2)
  map2(n1, n2, ~ extract_shortest_path(gr, .x, .y))
}


node_pal <- c(
  "source" = "dodgerblue",
  "destination" = "tomato",
  "on path" = "grey30"
)

modify_shortest_graph <- function(gr, n1, n2) {
  gr %N>%
    mutate(
      node_color = case_when(
        name == !!n1 ~ "source",
        name == !!n2 ~ "destination",
        TRUE ~ "on path"
      ),
      node_color = factor(node_color, levels = names(node_pal))
    )
}


plot_shortest_path_gr <- function(gr, pal = node_pal) {
  p <- gr %>%
    ggraph("stress") +
    geom_edge_link(
      color = "grey40",
      alpha = 0.3,
      width = 0.4
    ) +
    geom_node_point(
      aes(color = node_color, size = node_color)
    ) +
    scale_color_manual(values = pal) +
    scale_size_manual(values = c(2, 2, 1)) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3,
      family = "Arial"
    ) +
    theme_graph() +
    theme(
      legend.title = element_blank()
    )
  return(p)
}


plot_shortest_paths_grs <- function(gr, n1, n2) {
  n1 <- match_length(n1, n2)
  for (i in 1:length(n1)) {
    mod_gr <- gr[[i]] %>%
      modify_shortest_graph(n1[[i]], n2[[i]]) %>%
      plot_shortest_path_gr()

    ggsave_wrapper(
      mod_gr,
      plot_path(
        GRAPHS_DIR,
        glue("shortest-path_{n1[[i]]}-{n2[[i]]}.svg")
      ),
      "medium"
    )
  }
  return(NULL)
}


pretty_print_results <- function(tib) {
  col_names <- c(
    "Gene",
    "Shortest dist.",
    "Avg. shortest dist.", "std. dev", "p-value",
    "Num. of shortest paths",
    "Avg. number paths", "std. dev.", "p-value"
  )
  tib %>%
    select(
      hugo_symbol,
      len_sp, avg_len_sp, sd_len_sp, len_pval,
      num_sp, avg_num_sp, sd_num_sp, num_pval
    ) %>%
    mutate(
      sd_len_sp = round(sd_len_sp, 1),
      sd_num_sp = round(sd_num_sp, 1),
      len_pval = round(len_pval, 3),
      num_pval = round(num_pval, 3)
    ) %>%
    data.table::setnames(col_names) %T>%
    write_tsv(table_path(GRAPHS_DIR, "shortest-path-statistics.tsv")) %>%
    knitr::kable(col.names = col_names) %>%
    print()
  return(tib)
}


# `side = "greater"`: p-value for x is larger than random
# `side = "less"`: p-value for x is smaller than random
calculate_bootstrap_pvalues <- function(x, rs, side = c("greater", "less")) {
  rs <- map(rs, ~ .x[is.finite(.x)])
  pvals <- NULL
  if (side == "greater") {
    pvals <- map2_dbl(x, rs, ~ 1 - sum(.x > .y) / length(.y))
  } else if (side == "less") {
    pvals <- map2_dbl(x, rs, ~ 1 - sum(.x < .y) / length(.y))
  } else {
    stop(glue("Size '{side}' is not recognized."))
  }
  return(pvals)
}


tibble(hugo_symbol = YJ_GENES) %>%
  mutate(
    len_sp = shortest_distances(combined_ppi_gr, "KRAS", hugo_symbol),
    num_sp = number_shortest_paths(combined_ppi_gr, "KRAS", hugo_symbol),
    sp_gr = extract_shortest_paths(combined_ppi_gr, "KRAS", hugo_symbol),
    sp_plt = plot_shortest_paths_grs(sp_gr, "KRAS", hugo_symbol),
    rdm_sp = random_shortest_distances(combined_ppi_gr, hugo_symbol),
    avg_len_sp = map_dbl(rdm_sp, mean),
    sd_len_sp = map_dbl(rdm_sp, sd),
    len_pval = calculate_bootstrap_pvalues(len_sp, rdm_sp, "less"),
    rdm_num_sp = random_number_shortest_paths(combined_ppi_gr, hugo_symbol),
    avg_num_sp = map_dbl(rdm_num_sp, mean),
    sd_num_sp = map_dbl(rdm_num_sp, sd),
    num_pval = calculate_bootstrap_pvalues(len_sp, rdm_num_sp, "greater"),
  ) %>%
  pretty_print_results()


#### ---- Follow-up ---- ####

combined_ppi_gr %N>%
  filter(name %in% c("KRAS", "BAG4", "HSPA13"))

string_gr %N>%
  filter(name %in% c("KRAS", "BAG4", "HSPA13")) %E>%
  pull(.orig_data) %>%
  bind_rows() %>%
  glimpse()



#### ---- G13D comutation for Yi-Jang ---- ####

out_file <- "COAD_G13D_comutation-genes_all-mutations_annotated.tsv"

coad_g13d_comut_genes <- genetic_interaction_df %>%
  filter(kras_allele == "KRAS_G13D" & cancer == "COAD") %>%
  select(hugo_symbol, genetic_interaction, interaction_p_value = p_val)

coad_g13d_comut_annotated <- cancer_coding_av_muts_df %>%
  select(-ras) %>%
  filter(cancer == "COAD") %>%
  right_join(coad_g13d_comut_genes, by = "hugo_symbol")

write_tsv(
  coad_g13d_comut_annotated,
  table_path(GRAPHS_DIR, out_file)
)
