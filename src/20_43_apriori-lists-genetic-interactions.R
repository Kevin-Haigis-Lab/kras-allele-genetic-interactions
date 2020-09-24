
# Checking for genes-of-interest (goi) in the genetic interactions with KRAS

GRAPHS_DIR <- "20_43_apriori-lists-genetic-interactions"
reset_graph_directory(GRAPHS_DIR)


#### ---- FDR analysis of comuts. in a priori gene sets ---- ####

ProjectTemplate::cache("wide_genetic_interaction_df",
  depends = "genetic_interaction_gr",
  {
    wide_genetic_interaction_df <- kegg_geneset_df %>%
      group_by(hugo_symbol) %>%
      summarise(KEGG = paste(gene_set, collapse = ", ")) %>%
      ungroup() %>%
      full_join(
        cosmic_cgc_df %>%
          distinct(hugo_symbol) %>%
          add_column(CGC = TRUE),
        by = "hugo_symbol"
      ) %>%
      full_join(
        kras_interactors_bioid_df %>%
          distinct(hugo_symbol) %>%
          add_column(BioID = TRUE),
        by = "hugo_symbol"
      )
    return(wide_genetic_interaction_df)
  }
)


apriori_genes_fdr_adjusted_comutations <- complete_genetic_interaction_df %>%
  inner_join(
    wide_genetic_interaction_df,
    by = "hugo_symbol"
  ) %>%
  mutate(
    KEGG = !is.na(KEGG),
    CGC = !is.na(CGC),
    BioID = !is.na(BioID)
  ) %>%
  select(hugo_symbol:p_val, genetic_interaction, is_sig, KEGG:BioID) %>%
  pivot_longer(
    c(KEGG, CGC, BioID),
    names_to = "gene_set",
    values_to = "is_in_gene_set"
  ) %>%
  filter(is_in_gene_set) %>%
  select(-is_in_gene_set) %>%
  group_by(cancer, kras_allele, gene_set, genetic_interaction) %>%
  mutate(bh_adj_p_val = p.adjust(p_val, method = "BH")) %>%
  ungroup()



#### ---- Plot comutation networks ---- ####


# Queryable list of which plots to save to protos for Figure 2.
imgs_to_save_for_figure <- list(
  COAD = list(suffix = "_allLists"),
  LUAD = list(suffix = c("_allLists", "_kegg")),
  PAAD = list(suffix = "_allLists")
)


filter_for_adjusted_p_values <- function(gr,
                                         gene_sets = "all",
                                         adj_val_cutoff = 0.25) {
  nodes <- as_tibble(gr, active = "nodes")

  cancer <- as_tibble(gr, active = "edges") %>%
    pull(cancer) %>%
    unique()

  sig_interactions <- apriori_genes_fdr_adjusted_comutations %>%
    filter(bh_adj_p_val < !!adj_val_cutoff & is_sig) %>%
    filter(cancer == !!cancer) %>%
    filter(gene_set %in% !!gene_sets | !!gene_sets == "all") %>%
    distinct(hugo_symbol, kras_allele, genetic_interaction)

  gr %E>%
    inner_join(
      sig_interactions,
      by = c("hugo_symbol", "kras_allele", "genetic_interaction")
    )
}



# Manual adjustments to some nodes in the plots for figures.
adjust_layout_manually <- function(layout, CANCER, SUFFIX) {
  layout_attrs <- attributes(layout)

  msg <- glue("Manually adjust plot for {CANCER} with suffix '{SUFFIX}'")


  if (CANCER == "COAD" & SUFFIX == "_allLists") {
    message(msg)

    shift_down <- c(
      "SMAD3", "G12S", "CACNA1E", "Q61H", "LRP2", "Q61L", "BRCA2", "TTN",
      "A146V", "MCM4", "G12C", "MTOR"
    )

    shift_left <- c(
      c("TTN", "A146V", "MCM4", "G12C", "MTOR")
    )

    layout <- layout %>%
      mutate(
        x = ifelse(node_label == "RICTOR", x[node_label == "G13D"], x),
        x = ifelse(node_label %in% shift_left, x - 0.6, x),
        y = ifelse(node_label %in% shift_down, y - 1.3, y)
      )
  }

  if (CANCER == "LUAD" & SUFFIX == "_kegg") {
    message(msg)
    layout <- layout %>%
      mutate(
        x = ifelse(node_label == "G13C", x + 1, x),
        x = ifelse(node_label == "PLCB3", x[node_label == "G13C"], x),
        y = ifelse(node_label == "PLCB3", y + 0.8, y),
        x = ifelse(node_label == "EGFR", x - 0.7, x),
        y = ifelse(node_label == "EGFR", y + 0.5, y),
        x = ifelse(node_label == "TP53", x - 0.3, x),
        y = ifelse(node_label == "TP53", y + 0.3, y)
      )
  }

  attributes(layout) <- layout_attrs
  return(layout)
}


plot_genetic_interaction_graph <- function(gr_to_plot, CANCER, SUFFIX = "",
                                           gr_layout = "fr") {
  set.seed(0)
  num_nodes <- igraph::vcount(gr_to_plot)

  layout <- create_layout(
    gr_to_plot,
    layout = "fr",
    start.temp = igraph::vcount(gr_to_plot)^(1 / 4)
  )

  layout <- adjust_layout_manually(layout, CANCER, SUFFIX)

  gr_plot <- ggraph(layout) +
    geom_edge_link(
      aes(
        color = genetic_interaction,
        width = -log(p_val + 0.0000001)
      )
    ) +
    scale_edge_color_manual(
      values = comut_mutex_pal,
      guide = FALSE
    ) +
    scale_edge_width_continuous(
      range = c(0.2, 1.5),
      guide = guide_legend(
        label.position = "top",
        keyheight = unit(1, "mm")
      )
    ) +
    geom_node_label(
      aes(
        label = node_label,
        fontface = label_face,
        fill = node_fill,
        color = node_color,
        size = node_size
      ),
      family = "Arial",
      label.padding = unit(0.07, "lines"),
      label.r = unit(0.1, "lines"),
      label.size = 0
    ) +
    scale_color_identity() +
    scale_fill_manual(
      values = short_allele_pal,
      guide = FALSE,
      na.value = "grey85"
    ) +
    scale_size_manual(
      values = c(big = 1.6, small = 1.2),
      guide = FALSE
    ) +
    theme_graph() +
    theme(
      text = element_text(family = "Arial")
    ) +
    labs(
      edge_width = "-log( p-value )"
    )
  save_path <- plot_path(
    GRAPHS_DIR,
    glue("goi_overlap_genetic_interactions_network_{CANCER}{SUFFIX}.svg")
  )
  ggsave_wrapper(gr_plot, save_path, "small")

  fig_info <- imgs_to_save_for_figure[[CANCER]]
  if (!is.null(fig_info) & SUFFIX %in% fig_info$suffix) {
    base_n <- file_sans_ext(basename(save_path))
    saveFigRds(gr_plot, base_n)
  }
}


for (CANCER in unique(genetic_interaction_df$cancer)) {
  gr_to_plot <- genetic_interaction_gr %N>%
    left_join(wide_genetic_interaction_df,
      by = c("name" = "hugo_symbol")
    ) %E>%
    filter(cancer == !!CANCER) %N>%
    filter(is_kras | !is.na(KEGG) | !is.na(CGC) | !is.na(BioID)) %>%
    filter_for_adjusted_p_values() %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(
      label_face = ifelse(is_kras, "bold", "plain"),
      node_label = str_remove_all(name, "KRAS_"),
      node_fill = ifelse(is_kras, node_label, NA),
      node_color = ifelse(node_label %in% kras_dark_lbls,
        "white", "black"
      ),
      node_size = ifelse(is_kras, "big", "small")
    )

  if (igraph::vcount(gr_to_plot) == 0) {
    next
  }

  # plot all interactions with goi
  plot_genetic_interaction_graph(gr_to_plot, CANCER, "_allLists")

  plot_specific_lists <- function(gene_list, suffix) {
    j <- which(colnames(as_tibble(gr_to_plot, "nodes")) == gene_list)
    idx <- !is.na(as_tibble(gr_to_plot, "nodes")[, j])

    gr_to_plot_MOD <- gr_to_plot %N>%
      filter(is_kras | !!idx) %>%
      filter(centrality_degree(mode = "all") > 0)
    if (igraph::vcount(gr_to_plot_MOD) > 0) {
      grlay <- ifelse(CANCER == "LUAD" & suffix == "_kegg",
        "stress", "fr"
      )

      plot_genetic_interaction_graph(gr_to_plot_MOD, CANCER, suffix,
        gr_layout = grlay
      )
    }
  }

  tibble(
    gene_list = c("KEGG", "CGC", "BioID"),
    suffix = c("_kegg", "_cgc", "_BioID")
  ) %>%
    pmap(plot_specific_lists)
}
