# Analyzing the results of GSEA on the DepMap data.

GRAPHS_DIR <- "10_37_gsea-depmap-analysis"
reset_graph_directory(GRAPHS_DIR)


#### ---- Read in results ---- ####

gsea_df <- file.path("data", "gsea", "output") %>%
  list.dirs(recursive = FALSE) %>%
  tibble(dir = .) %>%
  mutate(
    dir_base = basename(dir),
    cancer = str_extract(dir_base, "^[:alpha:]+(?=_)"),
    allele = str_extract(dir_base, "(?<=_)[:alnum:]+(?=\\.G)"),
    timestamp = str_extract(dir_base, "(?<=Gsea\\.)[:digit:]+$"),
    timestamp = lubridate::as_datetime(as.numeric(timestamp) / 1000),
    data = get_gsea_reports(dir)
  ) %>%
  unnest(data) %>%
  mutate(
    gene_set_family = str_split_fixed(name, "_", 2)[, 1],
    gene_set = str_split_fixed(name, "_", 2)[, 2]
  ) %>%
  mutate(
    nes = ifelse(cancer == "PAAD" & allele == "G12V", -nes, nes),
    es = ifelse(cancer == "PAAD" & allele == "G12V", -es, es)
  )

ProjectTemplate::cache("gsea_df")


# Pull an enrichment plot ("enplot") from the GSEA results for a the specified
#   cancer, allele, and name (gene-set).
ENPLOT_GRAPHS_DIR <- "10_37_gsea-depmap-output"
reset_graph_directory(ENPLOT_GRAPHS_DIR)

pull_gsea_enplot <- function(cancer, allele, gene_set, gene_set_family,
                             force_overwrite = FALSE,
                             ...) {
  cancer_allele <- paste0(cancer, "_", allele)

  gsea_dirs <- file.path("data", "gsea", "output") %>%
    list.dirs(recursive = FALSE)
  gsea_dir <- gsea_dirs[str_detect(gsea_dirs, cancer_allele)]

  svg_files <- list.files(gsea_dir, pattern = "^enplot.*svg.gz", full.names = TRUE)

  if (str_length(gene_set) > 100) {
    gene_set_regex <- glue("{gene_set_family}_{str_sub(gene_set, 1, 100)}")
  } else {
    gene_set_regex <- glue("{gene_set_family}_{gene_set}_[:digit:]+")
  }

  svg_file <- svg_files[str_detect(basename(svg_files), gene_set_regex)]

  target_dir <- plot_path(ENPLOT_GRAPHS_DIR, cancer_allele)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir)
  }

  target_path <- file.path(target_dir, basename(svg_file))
  unziped_target_path <- gsub("[.]gz$", "", target_path)

  if (length(target_path) != 1) {
    return(NULL)
  }

  if (file.exists(unziped_target_path) & !force_overwrite) {
    return(NULL)
  }

  a <- file.copy(svg_file, target_path)
  a <- R.utils::gunzip(target_path, overwrite = TRUE)

  return(NULL)
}


#### ---- Plot GSEA results ---- ####


uninteresting_terms_regex <- c(
  "pancreas", "keratin", "disease", "muscle"
) %>%
  paste0(collapse = "|") %>%
  regex(ignore_case = TRUE)


standardize_names <- function(x) {
  str_to_sentence(x) %>%
    str_replace_all("_", " ")
}


gsea_plot <- function(tib, title_suffix = "") {
  p <- tib %>%
    ggplot() +
    geom_point(
      aes(
        x = allele,
        y = gene_set,
        color = nes,
        size = -log10(fdr_q_val)
      )
    ) +
    scale_color_gradient2(
      low = synthetic_lethal_pal["down"],
      high = synthetic_lethal_pal["up"]
    ) +
    theme_bw() +
    theme(
      text = element_text("arial"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      title = glue("GSEA of Alleles in {title_suffix}"),
      color = "NES",
      size = "-log10( adj. p-val. )"
    )
  return(p)
}


plot_gsea_results <- function(cancer, data) {
  mod_data <- standard_gsea_results_filter(data)

  if (nrow(mod_data) == 0) {
    return()
  }

  mod_data <- mod_data %>%
    mutate(
      gene_set = standardize_names(gene_set),
      gene_set = str_wrap(gene_set, width = 30)
    )

  p <- mod_data %>%
    gsea_plot(title_suffix = cancer)
  save_path <- plot_path(GRAPHS_DIR, glue("gsea-results-{cancer}-all.svg"))
  ggsave_wrapper(p, save_path, "large")

  for (fam in unique(data$gene_set_family)) {
    suffix <- glue("{cancer} ({fam})")
    tt <- mod_data %>%
      filter(gene_set_family == !!fam)

    if (nrow(tt) == 0) {
      next
    }

    p <- gsea_plot(tt, title_suffix = suffix)
    save_path <- plot_path(
      GRAPHS_DIR,
      glue("gsea-results-{cancer}-{fam}.svg")
    )
    ggsave_wrapper(p, save_path, "large")
  }
}


gsea_df %>%
  filter(gene_set_family %in% c(
    "HALLMARK", "KEGG", "REACTOME",
    "BIOCARTA", "PID"
  )) %>%
  group_by(cancer) %>%
  nest() %>%
  purrr::pwalk(plot_gsea_results)



#### ---- Plot select GSEA results ---- ####

select_gsea_results <- tibble::tribble(
  ~cancer, ~gene_set_family, ~gene_set, ~gene_set_pretty,
  "COAD", "REACTOME", "Tp53 regulates metabolic genes", "TP53 regulates metabolic genes",
  "COAD", "REACTOME", "Respiratory electron transport", "Respiratory electron transport",
  "COAD", "REACTOME", "Regulation of expression of slits and robos",
  "Regulation of expression of SLITS and ROBOS",
  "COAD", "REACTOME", "Nonsense mediated decay nmd", "Nonsense Mediated Decay",
  "COAD", "REACTOME", "Complex i biogenesis", "Complex I biogenesis",
  "COAD", "REACTOME", "Mitotic spindle checkpoint", "Mitotic spindle checkpoint",
  "COAD", "REACTOME", "Downstream signaling of activated fgfr1",
  "Downstream signaling of activated FGFR1",
  "COAD", "REACTOME", "Resolution of d loop structures through synthesis dependent strand annealing sdsa",
  "Resolution of D-loop structures through SDSA",
  "COAD", "REACTOME", "Fanconi anemia pathway", "Fanconi anemia pathway",
  "COAD", "REACTOME", "Dual incision in tc ner", "Dual incision in TC-NER",
  "COAD", "REACTOME", "Tgf beta receptor signaling in emt epithelial to mesenchymal transition",
  "TGF-β receptor signaling in EMT",
  "COAD", "REACTOME", "Cell extracellular matrix interactions", "ECM interactions",
  "COAD", "HALLMARK", "Oxidative phosphorylation", "Oxidative phosphorylation",
  "COAD", "BIOCARTA", "Tff pathway", "Trefoil factor (TFF) pathway",
  "COAD", "KEGG", "Mismatch repair", "Mismatch repair",
  "COAD", "PID", "Erbb4 pathway", "ERBB4 pathway",
  "COAD", "PID", "Bard1 pathway", "BARD1 pathway",

  "PAAD", "REACTOME", "G2 m dna damage checkpoint", "G2/M DNA damage checkpoint",
  "PAAD", "REACTOME", "Mitochondrial translation", "Mitochondrial translation",
  "PAAD", "REACTOME", "Activation of gene expression by srebf srebp",
  "Activation of gene expression by SREBF and SREBP",
  "PAAD", "REACTOME", "Tp53 regulates metabolic genes", "TP53 regulates metabolic genes",
  "PAAD", "BIOCARTA", "Toll pathway", "Toll pathway",
  "PAAD", "PID", "Fanconi pathway", "Fanconi pathway",
  "PAAD", "PID", "Foxm1 pathway", "FOXM1 pathway",
  "PAAD", "PID", "Atm pathway", "ATM pathway",
  "PAAD", "PID", "S1p s1p2 pathway", "S1P S1P2 pathway",
  "PAAD", "KEGG", "Peroxisome", "Peroxisome",
  "PAAD", "REACTOME", "Pi 3k cascade:fgfr1", "PI-3K cascade:FGFR1",
  "PAAD", "REACTOME", "Nonhomologous end joining nhej", "Nonhomologous end joining",
  "PAAD", "REACTOME", "Homology directed repair", "Homology directed repair",
  "PAAD", "REACTOME", "Fgfr1 mutant receptor activation", "FGFR1 mutant receptor activation",
  "PAAD", "REACTOME", "Fanconi anemia pathway", "Fanconi anemia pathway",
  "PAAD", "REACTOME", "Cellular senescence", "Cellular senescence",
  "PAAD", "REACTOME", "Signaling by egfr in cancer", "Signaling by EGFR in cancer",
  "PAAD", "REACTOME", "Wnt5a dependent internalization of fzd4", "WNT5A-dependent internalization of FZD4",
  "PAAD", "REACTOME", "Resolution of d loop structures", "Resolution of D-loop structures"
)


filter_gsea_for_select_results <- function(cancer, df, key) {
  key <- key %>%
    filter(cancer == !!cancer)
  df %>%
    mutate(gene_set = standardize_names(gene_set)) %>%
    inner_join(key, by = c("gene_set_family", "gene_set")) %>%
    mutate(gene_set = gene_set_pretty)
}


# Information for where to save the ggplot proto object for a Figure.
cancers_to_save_for_figures <- c("COAD", "PAAD")


# MEMO sort of gene sets based on FDR and x-axis order.
get_geneset_order <- function(gene_set, x_order, y_metric) {
  geneset_order <- bind_cols(
    gene_set = gene_set,
    x_order = x_order,
    y_metric = y_metric
  ) %>%
    mutate(score = x_order**2 * y_metric) %>%
    group_by(gene_set) %>%
    summarize(score = sum(score)) %>%
    ungroup() %>%
    arrange(-score) %>%
    pull(gene_set)
  return(factor(gene_set, levels = geneset_order))
}


# Make dot-plots of selected gene sets.
select_gsea_plot <- function(cancer, data, ...) {
  mod_data <- standard_gsea_results_filter(data)

  if (nrow(mod_data) == 0) {
    return()
  }
  mod_data <- mod_data %>%
    mutate(
      allele = factor_alleles(allele),
      gene_set = str_wrap(gene_set, width = 40),
      gene_set = get_geneset_order(
        gene_set,
        as.numeric(allele),
        -log10(fdr_q_val)
      )
    )

  p <- gsea_plot(mod_data, title_suffix = cancer)
  save_path <- plot_path(
    GRAPHS_DIR,
    glue("gsea-results-{cancer}-select.svg")
  )
  ggsave_wrapper(p, save_path, "wide")

  purrr::pwalk(mod_data, pull_gsea_enplot)

  if (cancer %in% cancers_to_save_for_figures) {
    saveFigRds(p, basename(save_path))
  }
}


gsea_df %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(data = purrr::map2(cancer, data,
    filter_gsea_for_select_results,
    key = select_gsea_results
  )) %>%
  purrr::pwalk(select_gsea_plot)


#### ---- Ranking heatmap ---- ####

# Gene sets used in the GSEA
gsea_geneset_df <- bind_rows(msigdb_hallmark_df, msigdb_c2_df)


gsea_input_data <- depmap_modelling_df %>%
  filter_depmap_by_allele_count() %>%
  group_by(cancer) %>%
  filter(n_distinct(kras_allele) >= 3) %>%
  filter(!is_deleted) %>%
  add_count(cancer, hugo_symbol, kras_allele) %>%
  group_by(cancer, hugo_symbol) %>%
  filter(all(n >= 3)) %>%
  select(-n) %>%
  ungroup()

# Remove genes that do not have gene effect data for each cell line
gsea_input_data %<>%
  group_by(cancer, hugo_symbol) %>%
  mutate(num_cell_lines = n_distinct(dep_map_id)) %>%
  group_by(cancer) %>%
  filter(num_cell_lines == max(num_cell_lines)) %>%
  ungroup()


rank_depmap_data <- function(data) {
  data %>%
    group_by(hugo_symbol) %>%
    mutate(effect_rank = rank(gene_effect, ties.method = "random")) %>%
    ungroup()
}


get_alpha_values_by_distance <- function(data) {
  data %>%
    group_by(hugo_symbol) %>%
    mutate(
      alpha_val = abs(effect_rank - median(effect_rank)),
      alpha_val = scales::rescale(alpha_val, to = c(0.3, 1.0))
    ) %>%
    ungroup()
}


get_alpha_values_to_highlight_allele <- function(data, allele) {
  data %>%
    mutate(alpha_val = ifelse(kras_allele == !!allele, 1.0, 0.7))
}


get_color_values_to_highlight_allele <- function(data, allele) {
  data %>%
    mutate(color_val = ifelse(kras_allele == !!allele, "black", NA))
}


save_to_proto <- function(cancer, allele, geneset, gg_obj, save_name) {
  if (!(cancer %in% cancers_to_save_for_figures)) {
    return(NULL)
  }

  selected_genesets <- select_gsea_results %>%
    filter(cancer == !!cancer) %>%
    pull(gene_set) %>%
    standardize_names() %>%
    str_to_lower()
  cond <- str_detect(
    str_to_lower(standardize_names(geneset)),
    selected_genesets
  )
  if (any(cond)) {
    saveFigRds(gg_obj, basename(save_name))
  }
}


plot_ranked_data <- function(df, cancer, allele, geneset,
                             use_alpha_grad = FALSE) {
  plot_title <- str_replace_all(geneset, "_", " ") %>%
    str_to_sentence() %>%
    str_wrap(50)

  p <- df %>%
    get_alpha_values_by_distance() %>%
    get_color_values_to_highlight_allele(allele = allele) %>%
    ggplot(aes(x = effect_rank, y = hugo_symbol))

  if (use_alpha_grad) {
    p <- p +
      geom_tile(aes(fill = kras_allele, alpha = alpha_val), size = 0.5)
  } else {
    p <- p +
      geom_tile(aes(fill = kras_allele), size = 0.5)
  }
  p <- p +
    scale_fill_manual(values = short_allele_pal) +
    scale_alpha_identity() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      axis.title = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.key.size = unit(2, "mm"),
      axis.text.y = element_text(face = "italic")
    ) +
    labs(
      title = glue("{cancer} - {allele}\n{plot_title}")
    )

  save_name <- plot_path(
    GRAPHS_DIR,
    glue("rankplot_{cancer}_{allele}_{geneset}.svg")
  )
  ggsave_wrapper(p, save_name, width = 5, height = 3)
  save_to_proto(cancer, allele, geneset, p, save_name)

  return(df)
}


get_genes_to_plot <- function(gsea_xls_data, cancer, n_genes) {
  gtp <- gsea_xls_data %>%
    slice(seq(1, n_genes)) %>%
    pull(probe) %>%
    rev()
}


plot_enrichment_heatmap <- function(cancer, name, allele, n_genes = 10, ...) {
  genes_to_plot <- get_geneset_enrichment_results(cancer, allele, name)

  if (is.null(genes_to_plot)) {
    return(NULL)
  }

  genes_to_plot <- get_genes_to_plot(genes_to_plot, cancer, n_genes)
  gsea_input_data %>%
    filter(cancer == !!cancer & hugo_symbol %in% !!genes_to_plot) %>%
    mutate(hugo_symbol = factor(hugo_symbol, levels = genes_to_plot)) %>%
    rank_depmap_data() %>%
    plot_ranked_data(
      cancer = cancer,
      allele = allele,
      geneset = name
    )
}


# Calculate a measure of enrichment for each allele at each rank.
calculate_enrichment_average <- function(df, rescale = TRUE) {
  n_rows <- n_distinct(df$hugo_symbol)
  df %>%
    group_by(kras_allele) %>%
    mutate(num_celllines = n_distinct(dep_map_id)) %>%
    ungroup() %>%
    group_by(effect_rank, kras_allele) %>%
    mutate(
      fraction_rows = n_distinct(hugo_symbol) / !!n_rows / num_celllines
    ) %>%
    group_by(effect_rank) %>%
    mutate(
      fraction_rows = ifelse(
        !!rescale,
        scales::rescale(fraction_rows, to = c(0.1, 0.8)),
        fraction_rows
      )
    ) %>%
    ungroup()
}


plot_ranked_bar <- function(df, cancer, allele, geneset) {
  plot_title <- str_replace_all(geneset, "_", " ") %>%
    str_to_sentence() %>%
    str_wrap(50)

  p <- df %>%
    ggplot(aes(x = effect_rank, y = kras_allele)) +
    geom_tile(aes(fill = kras_allele, alpha = fraction_rows), size = 0.5) +
    scale_fill_manual(values = short_allele_pal) +
    scale_alpha_identity() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      axis.title = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.key.size = unit(2, "mm")
    ) +
    labs(
      title = glue("{cancer} - {allele}\n{plot_title}")
    )

  save_name <- plot_path(
    GRAPHS_DIR,
    glue("rankbar_{cancer}_{allele}_{geneset}.svg")
  )
  ggsave_wrapper(p, save_name, width = 5, height = 3)
  save_to_proto(cancer, allele, geneset, p, save_name)

  return(df)
}


plot_ranked_density <- function(df, cancer, allele, geneset, n_ranks) {
  plot_title <- str_replace_all(geneset, "_", " ") %>%
    str_to_sentence() %>%
    str_wrap(50)

  p <- df %>%
    ggplot(aes(x = effect_rank)) +
    geom_density(aes(color = kras_allele, fill = kras_allele),
      size = 0.5, alpha = 0.2
    ) +
    scale_color_manual(values = short_allele_pal) +
    scale_fill_manual(values = short_allele_pal) +
    scale_x_continuous(
      limits = c(0, n_ranks),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      axis.title = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.key.size = unit(2, "mm")
    ) +
    labs(
      title = glue("{cancer} - {allele}\n{plot_title}")
    )

  save_name <- plot_path(
    GRAPHS_DIR, glue("rankline_{cancer}_{allele}_{geneset}.svg")
  )
  ggsave_wrapper(p, save_name, width = 5, height = 3)
  save_to_proto(cancer, allele, geneset, p, save_name)

  return(df)
}


plot_enrichment_bar <- function(cancer, name, allele, n_genes = 10, ...) {
  genes_to_plot <- get_geneset_enrichment_results(cancer, allele, name)

  if (is.null(genes_to_plot)) {
    return(NULL)
  }

  genes_to_plot <- get_genes_to_plot(genes_to_plot, cancer, n_genes)

  gsea_input_data %>%
    filter(cancer == !!cancer & hugo_symbol %in% !!genes_to_plot) %>%
    mutate(hugo_symbol = factor(hugo_symbol, levels = genes_to_plot)) %>%
    rank_depmap_data() %>%
    calculate_enrichment_average(rescale = TRUE) %>%
    mutate(kras_allele = fct_rev(kras_allele)) %T>%
    plot_ranked_bar(
      cancer = cancer,
      allele = allele,
      geneset = name
    )

  return(NULL)
}


plot_enrichment_density <- function(cancer, name, allele, n_genes = 10, ...) {
  genes_to_plot <- get_geneset_enrichment_results(cancer, allele, name)

  if (is.null(genes_to_plot)) {
    return(NULL)
  }

  genes_to_plot <- get_genes_to_plot(genes_to_plot, cancer, n_genes)

  dat <- gsea_input_data %>%
    filter(cancer == !!cancer & hugo_symbol %in% !!genes_to_plot) %>%
    mutate(hugo_symbol = factor(hugo_symbol, levels = genes_to_plot))
  dat %>%
    rank_depmap_data() %>%
    plot_ranked_density(
      cancer = cancer,
      allele = allele,
      geneset = name,
      n_ranks = n_distinct(dat$dep_map_id)
    )
}


gsea_df %>%
  standard_gsea_results_filter() %T>%
  pwalk(plot_enrichment_heatmap) %T>%
  pwalk(plot_enrichment_bar) %T>%
  pwalk(plot_enrichment_density)
