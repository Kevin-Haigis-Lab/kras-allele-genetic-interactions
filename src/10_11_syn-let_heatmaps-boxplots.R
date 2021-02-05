# Box-plots and heatmaps of results from linear modeling of allele-specific
# synthetic lethality.

set.seed(0)

#### ---- Box-plots ---- ####

# directory for save box-plots
GRAPHS_DIR_BOXES <- "10_11_linear-modeling-syn-let_boxplots"
reset_graph_directory(GRAPHS_DIR_BOXES)


save_boxplot_proto <- function(gg_obj, save_path, cancer) {
  saveFigRds(gg_obj, basename(save_path))
}


# plot the results of the first analysis
plot_pairwise_test_results <- function(hugo_symbol, cancer, data, ...) {
  set.seed(0)

  data <- distinct(data)

  stat_tib <- compare_means(
    gene_effect ~ kras_allele,
    data = data,
    method = "wilcox.test", p.adjust.method = "BH"
  ) %>%
    filter(p.adj < 0.05)

  stat_bar_height <- 0.08
  stat_bar_y_positions <- c(max(data$gene_effect) + stat_bar_height)
  for (i in seq(1, nrow(stat_tib))) {
    stat_bar_y_positions <- c(
      stat_bar_y_positions,
      stat_bar_y_positions[(i - 1)] + stat_bar_height
    )
  }

  stat_tib$y.position <- stat_bar_y_positions

  p <- ggboxplot(
    data,
    x = "kras_allele",
    y = "gene_effect",
    color = "kras_allele",
    add = "jitter",
    size = 1
  ) +
    stat_pvalue_manual(stat_tib, label = "p.adj", family = "Arial") +
    scale_color_manual(values = short_allele_pal) +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    labs(y = "depletion effect")

  plot_fname <- plot_path(
    GRAPHS_DIR_BOXES,
    glue("{cancer}-{hugo_symbol}.svg")
  )
  ggsave_wrapper(p, plot_fname, size = "small")
}


# Plot the results of the first analysis.
# This version of the function uses my own box-plot creation function
#   located in "lib/stats-boxplot.R".
plot_pairwise_test_results2 <- function(hugo_symbol,
                                        cancer,
                                        data,
                                        ova_pairs,
                                        up_spacing = 0.02,
                                        dn_spacing = 0,
                                        save_proto = FALSE,
                                        replace_svg = TRUE,
                                        ...) {
  set.seed(0)

  data <- distinct(data) %>%
    mutate(
      x = fct_drop(factor_alleles(kras_allele)),
      y = gene_effect
    )

  stats_df <- ova_pairs %>%
    filter(adj_p_value < 0.05) %>%
    left_join(
      data %>%
         group_by(kras_allele) %>%
         filter(y == max(y)) %>%
         ungroup() %>%
         select(kras_allele, y, x),
      by = c("allele" = "kras_allele")
    ) %>%
    mutate(
      label = format_pvalue_label(adj_p_value, add_p = FALSE),
      label = paste0("p=<br>", label)
    )

  p <- stats_boxplot_boxplot(
    data,
    box_color = kras_allele,
    point_size = 0.25,
    point_alpha = 0.8,
    box_size = 0.4
  ) +
    ggtext::geom_richtext(
      aes(x = x, y = y, label = label),
      data = stats_df,
      size = 2,
      color = "black",
      hjust = 0.5,
      vjust = -0.3,
      family = "Arial",
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    scale_color_manual(values = short_allele_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    labs(y = "depletion effect")

  plot_fname <- plot_path(
    GRAPHS_DIR_BOXES,
    glue("{cancer}-{hugo_symbol}_extra.svg")
  )
  ggsave_wrapper(p, plot_fname, size = "small")

  if (save_proto) {
    save_boxplot_proto(p, plot_fname)
  }
}

# Select genes for figures.
select_gene_boxplots <- tibble::tribble(
  ~cancer, ~hugo_symbol, ~up_spacing, ~dn_spacing,
  "COAD", "FAF2", 0.06, 0.06,
  "COAD", "KNTC1", 0.04, 0.06,
  "COAD", "LIN7C", 0.01, 0.06,
  "COAD", "NKD1", 0.06, 0.06,
  "COAD", "STARD9", 0.01, 0.04,
  "COAD", "TFPT", 0.06, 0.06,
  "PAAD", "BRI3BP", 0.02, 0.06,
  "PAAD", "EGLN2", 0.03, 0.06,
  "PAAD", "JUN", 0.06, 0.06,
  "PAAD", "KHDRBS1", 0.04, 0.06,
  "PAAD", "MAPK8", 0.06, 0.04,
  "PAAD", "MAPKBP1", 0.06, 0.06,
  "PAAD", "SPC24", 0.06, 0.06,
  "PAAD", "ZZEF1", 0.06, 0.06,
  "PAAD", "ZNF701", 0.06, 0.06,
  "PAAD", "ZSCAN31", 0.06, 0.06
)

depmap_model_workflow_res %>%
  filter_depmap_model_workflow_res() %>%
  # pwalk(plot_pairwise_test_results) %>%
  inner_join(select_gene_boxplots, by = c("cancer", "hugo_symbol")) %>%
  pwalk(plot_pairwise_test_results2, save_proto = TRUE)


#### ---- Heatmaps ---- ####

# directory for save box-plots
GRAPHS_DIR_HEAT <- "10_11_linear-modeling-syn-let_pheatmaps"
reset_graph_directory(GRAPHS_DIR_HEAT)

cancer_pheatmap_manager <- list(
  COAD = list(
    col_cuts = 4,
    row_cuts = 4
  ),
  PAAD = list(
    col_cuts = 3,
    row_cuts = 3
  )
)


prep_pheatmap_df <- function(data, method = c("scale", "normalize")) {
  f <- function(x) {
    x
  }
  if (method == "scale") {
    f <- function(x) {
      scale(x)[, 1]
    }
  } else if (method == "normalize") {
    f <- function(x) {
      scales::rescale(x, to = c(-2, 2))
    }
  } else {
    stop(glue("method not a possible option: {method}"))
  }

  df <- data %>%
    group_by(hugo_symbol) %>%
    mutate(gene_effect_scaled = f(gene_effect)) %>%
    ungroup() %>%
    select(hugo_symbol, dep_map_id, gene_effect_scaled) %>%
    unique() %>%
    group_by(hugo_symbol, dep_map_id) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    pivot_wider(
      names_from = dep_map_id,
      values_from = gene_effect_scaled
    ) %>%
    column_to_rownames("hugo_symbol")
  invisible(df)
}



# Custom annotation legends using `ggplot2::geom_tile()`.
make_annotation_tile <- function(v, grp) {
  df <- tibble::enframe(v) %>% filter(!is.na(name))

  if (grp == "heat") {
    df %<>% mutate(name = as.numeric(name))
  } else {
    df %<>%
      mutate(name = factor(name, levels = rev(sort(unique(name)))))
  }

  p <- df %>%
    ggplot(aes(x = "", y = name)) +
    geom_tile(aes(fill = value), color = NA) +
    scale_fill_identity() +
    labs(y = grp)
  return(p)
}


# Save a proto for a figure.
save_pheatmap_proto <- function(cancer, ph, save_path,
                                anno_pal = NULL, heat_pal = NULL) {
  saveFigRds(ph, basename(save_path))

  anno_save_path <- function(pal_name) {
    paste0(file_sans_ext(basename(save_path)), "_", pal_name, ".svg")
  }

  if (!is.null(anno_pal)) {
    g_allele <- make_annotation_tile(anno_pal$kras_allele, "allele")
    p <- anno_save_path("allelepal")
    saveFigRds(g_allele, p)
    g_cls <- make_annotation_tile(anno_pal$cluster, "cluster")
    p <- anno_save_path("clusterpal")
    saveFigRds(g_cls, p)
  }

  if (!is.null(heat_pal)) {
    g_heat <- make_annotation_tile(heat_pal, "heat")
    p <- anno_save_path("heatpal")
    saveFigRds(g_heat, p)
  }
}


# Mappings from default cluster assignments to order shown in pheatmap.
cluster_number_map <- list(
  COAD = tibble::tribble(
    ~default_cluster, ~cluster,
    1, 2,
    2, 1,
    3, 3,
    4, 4
  ),
  PAAD = tibble::tribble(
    ~default_cluster, ~cluster,
    1, 3,
    2, 1,
    3, 2
  )
)


# Update the clusters of the final data frame to those manually assigned in
# `cluster_number_map`.
update_clusters <- function(tib) {
  mapping_df <- tibble::enframe(cluster_number_map, name = "cancer") %>%
    unnest(value)

  tib %>%
    dplyr::rename(default_cluster = gene_cls) %>%
    left_join(mapping_df, by = c("cancer", "default_cluster")) %>%
    select(-default_cluster) %>%
    dplyr::rename(gene_cls = cluster)
}


# Get color palette for clusters. It is derived from viridis.
cluster_color_pal <- function(n_vals) {
  cols <- viridis::viridis_pal()(n_vals)
  names(cols) <- seq(1, n_vals)
  return(cols)
}


apriori_genes <- unique(unlist(c(
  kegg_geneset_df$hugo_symbol,
  cosmic_cgc_df$hugo_symbol,
  kras_interactors_bioid_df$hugo_symbol
)))



# Custom row names of only a priori genesets.
labels_for_apriori_genes <- function(df) {
  ifelse(rownames(df) %in% apriori_genes, rownames(df), "")
}


# Plot pretty heatmaps for a cancer.
plot_cancer_heatmaps <- function(cancer, data, screen,
                                 merge_luad = TRUE,
                                 row_dist_method = "euclidean",
                                 col_dist_method = "euclidean",
                                 row_hclust_method = "complete",
                                 col_hclust_method = "complete",
                                 save_proto = TRUE) {
  set.seed(0)
  mod_data <- prep_pheatmap_df(data, "normalize")

  row_hclust <- hclust(dist(mod_data, method = row_dist_method),
    method = row_hclust_method
  )
  col_hclust <- hclust(dist(t(mod_data), method = col_dist_method),
    method = col_hclust_method
  )

  col_anno <- data %>%
    select(dep_map_id, kras_allele) %>%
    distinct() %>%
    rename(allele = kras_allele) %>%
    column_to_rownames("dep_map_id")

  row_anno <- cutree(
    row_hclust,
    k = cancer_pheatmap_manager[[cancer]]$row_cuts
  ) %>%
    enframe(value = "default_cluster") %>%
    left_join(cluster_number_map[[cancer]], by = "default_cluster") %>%
    mutate(cluster = factor(cluster, levels = sort(unique(cluster)))) %>%
    select(-default_cluster) %>%
    column_to_rownames("name")

  alleles <- as.character(sort(unique(col_anno$allele)))
  anno_pal <- list(
    allele = short_allele_pal[alleles],
    cluster = cluster_color_pal(max(as.numeric(row_anno$cluster)))
  )

  pal <- c(synthetic_lethal_pal["down"], "grey95", synthetic_lethal_pal["up"])
  pal <- colorRampPalette(pal)(7)

  ph <- pheatmap::pheatmap(
    mod_data,
    color = pal,
    cluster_rows = row_hclust,
    cluster_cols = col_hclust,
    annotation_col = col_anno,
    annotation_row = row_anno,
    annotation_colors = anno_pal,
    cutree_rows = cancer_pheatmap_manager[[cancer]]$row_cuts,
    cutree_cols = cancer_pheatmap_manager[[cancer]]$col_cuts,
    treeheight_row = 8,
    treeheight_col = 8,
    fontsize = 5,
    silent = TRUE,
    border_color = NA,
    fontfamily = "Arial",
    fontface = 1
  )

  ph <- italicize_pheatmap_rownames(ph)

  save_path <- plot_path(
    GRAPHS_DIR_HEAT,
    glue("{cancer}_{screen}_{row_dist_method}_{row_hclust_method}_pheatmap.svg")
  )
  save_pheatmap_svg(ph, save_path, width = 4, height = 6)

  names(pal) <- seq(min(mod_data, na.rm = TRUE),
    max(mod_data, na.rm = TRUE),
    length.out = length(pal)
  )
  if (save_proto) {
    save_pheatmap_proto(cancer, ph, save_path,
      anno_pal = anno_pal, heat_pal = pal
    )
  }
}


cluster_genes <- function(cancer, data,
                          row_dist_method = "euclidean",
                          row_hclust_method = "complete") {
  mod_data <- prep_pheatmap_df(data, method = "normalize")

  gene_hclust <- hclust(dist(mod_data, method = row_dist_method),
    method = row_hclust_method
  )
  gene_cls <- cutree(
    gene_hclust,
    k = cancer_pheatmap_manager[[cancer]]$row_cuts
  ) %>%
    enframe(name = "hugo_symbol", value = "gene_cls")

  return(gene_cls)
}


# Run `plot_cancer_heatmaps()` with a bunch of distance and clustering methods.
plot_cancer_heatmaps_multiple_methods <- function(cancer, data, screen,
                                                  merge_luad = TRUE,
                                                  methods_df) {
  pwalk(methods_df, plot_cancer_heatmaps,
    cancer = cancer, data = data,
    screen = screen, merge_luad = merge_luad
  )
}


# Combinations of `dist()` and `hclust()` methods.
dist_methods <- c("euclidean", "manhattan")
hclust_methods <- c("ward.D2", "single", "complete", "average")
methods_tib <- expand.grid(dist_methods, hclust_methods,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  set_names(c("row_dist_method", "row_hclust_method"))


# make a heatmap for the genes in each cancer
depmap_gene_clusters <- depmap_model_workflow_res %>%
  filter_depmap_model_workflow_res() %>%
  select(hugo_symbol, cancer, data) %>%
  unnest(data) %>%
  group_by(cancer) %>%
  nest() %T>%
  pwalk(plot_cancer_heatmaps_multiple_methods,
    screen = "CRISPR", methods_df = methods_tib
  ) %>%
  mutate(cluster_tib = purrr::map2(cancer, data, cluster_genes,
    row_dist_method = "manhattan",
    row_hclust_method = "ward.D2"
  )) %>%
  select(-data) %>%
  unnest(cluster_tib) %>%
  ungroup() %>%
  update_clusters()

cache("depmap_gene_clusters",
  depends = "depmap_model_workflow_res"
)


# Print out the number of genes per cancer.
depmap_gene_clusters %>%
  group_by(cancer) %>%
  summarise(num_genes = n_distinct(hugo_symbol)) %>%
  ungroup() %>%
  knitr::kable()
# > |cancer | num_genes|
# > |:------|---------:|
# > |COAD   |        62|
# > |PAAD   |       130|

#### ---- Supp. Data: heatmaps as numeric matrices ---- ####
# 3 decimal points

depmap_supp_data_df <- depmap_model_workflow_res %>%
  select(hugo_symbol, cancer, data) %>%
  right_join(depmap_gene_clusters, by = c("cancer", "hugo_symbol"))

supp_data_nums <- list(
  COAD = 8,
  PAAD = 9
)

for (cancer in names(supp_data_nums)) {
  depmap_supp_data_df %>%
    filter(cancer == !!cancer) %>%
    unnest(data) %>%
    mutate(gene_effect = scales::label_number(0.001)(gene_effect)) %>%
    mutate(cell_line_id = paste(kras_allele, dep_map_id, sep = " - ")) %>%
    select(cancer, hugo_symbol, gene_cls, cell_line_id, gene_effect) %>%
    pivot_wider(c(cancer, hugo_symbol, gene_cls),
      names_from = cell_line_id,
      values_from = gene_effect
    ) %>%
    arrange(gene_cls, hugo_symbol) %>%
    save_supp_data(
      supp_data_nums[[cancer]],
      glue("KRAS genetic dep in {cancer}")
    )
}
