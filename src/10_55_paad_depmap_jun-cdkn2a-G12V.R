# Investigating the relationship between JUN and KRAS G12V in PAAD.

GRAPHS_DIR <- "10_55_paad_depmap_jun-cdkn2a-G12V"
reset_graph_directory(GRAPHS_DIR)

jun_pw <- c(
  "KRAS", "TAB1", "MAP2K4", "MAP2K7", "MAPK8", "MAPK9", "MAPK10",
  "NR2C2", "MEN1", "JUN", "FOS", "ELK1", "JUND", "ATF2", "JDP2",
  "UBE2N", "USP28", "ZNF304", "CDT1", "KAT7",
  "CDKN2A", "CDKN2B", "TP53"
)


#### ---- Data preparation ---- ####

paad_cell_lines <- depmap_modelling_df %>%
  filter_depmap_by_allele_count() %>%
  group_by(cancer) %>%
  filter(n_distinct(kras_allele) >= 3) %>%
  ungroup() %>%
  filter(cancer == "PAAD") %>%
  select(dep_map_id, cancer, kras_allele) %>%
  unique()

jun_pw[!(jun_pw %in% depmap_modelling_df[depmap_modelling_df$cancer == "PAAD", ]$hugo_symbol)]


junpw_depmap <- gene_effect %>%
  inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
  filter(hugo_symbol %in% jun_pw)

all(jun_pw %in% junpw_depmap$hugo_symbol)
# jun_pw[!(jun_pw %in% junpw_depmap$hugo_symbol)]


junpw_cn <- ccle_copy_number %>%
  inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
  filter(hugo_symbol %in% jun_pw)

all(jun_pw %in% junpw_cn$hugo_symbol)


junpw_mut <- ccle_mutations %>%
  inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
  filter(hugo_symbol %in% jun_pw) %>%
  filter(hugo_symbol != "KRAS") %>%
  select(
    dep_map_id, hugo_symbol, variant_classification, variant_type,
    protein_change, is_deleterious, is_cosmic_hotspot,
    cancer, kras_allele
  ) %>%
  group_by(cancer, kras_allele, dep_map_id, hugo_symbol) %>%
  summarise(
    n_muts = n_distinct(protein_change),
    protein_change = paste(protein_change, collapse = ", "),
    is_deleterious = paste(is_deleterious, collapse = ", "),
    is_cosmic_hotspot = paste(is_cosmic_hotspot, collapse = ", ")
  ) %>%
  ungroup()

all(jun_pw %in% junpw_mut$hugo_symbol)


junpw_expr <- ccle_expression %>%
  inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
  filter(hugo_symbol %in% jun_pw)

all(jun_pw %in% junpw_expr$hugo_symbol)

cn_levels <- c("het_del", "norm", "amp")
joining_cols <- c("hugo_symbol", "dep_map_id", "cancer", "kras_allele")
junpw_depmap %<>%
  select(cancer, dep_map_id, kras_allele, hugo_symbol, gene_effect) %>%
  left_join(junpw_cn, by = joining_cols) %>%
  left_join(junpw_mut, by = joining_cols) %>%
  left_join(junpw_expr, by = joining_cols) %>%
  mutate(
    is_mutated = ifelse(!is.na(protein_change), "mut.", "not mut."),
    copy_number_label = factor(copy_number_label, levels = cn_levels)
  )

junpw_depmap %>%
  group_by(dep_map_id, hugo_symbol) %>%
  count() %>%
  filter(n > 1)


#### ---- Basic analysis ---- ####

set.seed(0)

# Box-plots of gene effect by KRAS allele
geneeffect_boxplots <- junpw_depmap %>%
  mutate(hugo_symbol = factor(hugo_symbol, levels = jun_pw)) %>%
  ggplot(aes(x = kras_allele, y = gene_effect)) +
  facet_wrap(. ~ hugo_symbol, scales = "free", ncol = 4) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey25", size = 0.6) +
  geom_boxplot(aes(color = kras_allele, fill = kras_allele),
    alpha = 0.5, outlier.shape = NA
  ) +
  geom_jitter(aes(color = kras_allele, shape = is_mutated),
    width = 0.25, size = 0.8
  ) +
  scale_color_manual(values = short_allele_pal) +
  scale_fill_manual(values = short_allele_pal) +
  scale_shape_manual(values = c(17, 16)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  geneeffect_boxplots,
  plot_path(GRAPHS_DIR, "geneeffect_boxplots.svg"),
  "large"
)


# Scatter plots of gene effect vs. RNA expr for each gene
geneeffect_rnaexpr_scatter <- junpw_depmap %>%
  mutate(hugo_symbol = factor(hugo_symbol, levels = jun_pw)) %>%
  ggplot(aes(x = rna_expression, y = gene_effect)) +
  facet_wrap(. ~ hugo_symbol, scales = "free", ncol = 4) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey50", size = 0.3) +
  geom_vline(xintercept = 0, linetype = 1, color = "grey50", size = 0.3) +
  geom_point(aes(
    color = kras_allele,
    shape = is_mutated,
    size = copy_number_label
  ),
  alpha = 0.6
  ) +
  scale_color_manual(values = short_allele_pal) +
  scale_shape_manual(values = c(17, 16)) +
  scale_size_manual(values = c(1, 2, 3)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  geneeffect_rnaexpr_scatter,
  plot_path(GRAPHS_DIR, "geneeffect_rnaexpr_scatter.svg"),
  "large"
)


# JUN expression and CDKN2A gene effect.
JUNexpr_CDKN2Adep_scatter <- junpw_depmap %>%
  filter(hugo_symbol %in% c("JUN", "CDKN2A")) %>%
  select(
    dep_map_id, kras_allele, hugo_symbol,
    gene_effect, rna_expression
  ) %>%
  pivot_wider(c(dep_map_id, kras_allele),
    names_from = hugo_symbol,
    values_from = c(rna_expression, gene_effect)
  ) %>%
  ggplot(aes(x = gene_effect_JUN, rna_expression_CDKN2A)) +
  geom_point(aes(color = kras_allele)) +
  scale_color_manual(values = short_allele_pal) +
  scale_shape_manual(values = c(17, 16)) +
  scale_size_manual(values = c(1, 2, 3)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  JUNexpr_CDKN2Adep_scatter,
  plot_path(GRAPHS_DIR, "JUNexpr_CDKN2Adep_scatter.svg"),
  "small"
)


# JUN gene effect vs. TP53 gene effect.
JUN_TP53_scatter <- junpw_depmap %>%
  filter(hugo_symbol %in% c("JUN", "TP53")) %>%
  select(
    dep_map_id, kras_allele, hugo_symbol,
    gene_effect, rna_expression
  ) %>%
  pivot_wider(c(dep_map_id, kras_allele),
    names_from = hugo_symbol,
    values_from = c(rna_expression, gene_effect)
  ) %>%
  ggplot(aes(x = gene_effect_JUN, gene_effect_TP53)) +
  geom_point(aes(color = kras_allele)) +
  scale_color_manual(values = short_allele_pal) +
  scale_shape_manual(values = c(17, 16)) +
  scale_size_manual(values = c(1, 2, 3)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  JUN_TP53_scatter,
  plot_path(GRAPHS_DIR, "JUN_TP53_scatter.svg"),
  "small"
)


# JUN gene effect vs. MEN1 gene effect.
JUN_MEN1_scatter <- junpw_depmap %>%
  filter(hugo_symbol %in% c("JUN", "MEN1")) %>%
  select(
    dep_map_id, kras_allele, hugo_symbol,
    gene_effect, rna_expression
  ) %>%
  pivot_wider(c(dep_map_id, kras_allele),
    names_from = hugo_symbol,
    values_from = c(rna_expression, gene_effect)
  ) %>%
  ggplot(aes(x = gene_effect_JUN, gene_effect_MEN1)) +
  geom_point(aes(color = kras_allele)) +
  scale_color_manual(values = short_allele_pal) +
  scale_shape_manual(values = c(17, 16)) +
  scale_size_manual(values = c(1, 2, 3)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  JUN_MEN1_scatter,
  plot_path(GRAPHS_DIR, "JUN_MEN1_scatter.svg"),
  "small"
)


#### ---- Comutation ---- ####

genetic_interaction_df %>%
  filter(cancer == "PAAD") %>%
  filter(hugo_symbol %in% jun_pw) %>%
  select(hugo_symbol, kras_allele, p_val, genetic_interaction)


#### ---- JUN-regulated genes ---- ####

jun_bs <- encode_tf_bindingsites %>%
  filter(gene_set == "JUN") %>%
  unique()

jun_bs

paad_deps <- depmap_model_workflow_res %>%
  filter_depmap_model_workflow_res() %>%
  filter(cancer == "PAAD")

jun_bs_deps <- paad_deps %>% filter(hugo_symbol %in% jun_bs$gene)

jun_bs_dep_boxplots <- jun_bs_deps %>%
  select(hugo_symbol, data) %>%
  unnest(data) %>%
  ggplot(aes(
    x = kras_allele, y = gene_effect,
    color = kras_allele, fill = kras_allele
  )) +
  facet_wrap(~hugo_symbol, scales = "free") +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_color_manual(values = short_allele_pal) +
  scale_fill_manual(values = short_allele_pal) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  jun_bs_dep_boxplots,
  plot_path(GRAPHS_DIR, "jun_bs_dep_boxplots.svg"),
  "large"
)



# Test for enrichment of JUN-regulated genes in the allele-specific
# genetic dependencies in PAAD.

all_paad_genes <- unique(depmap_modelling_df[depmap_modelling_df$cancer == "PAAD", ]$hugo_symbol)
sig_paad_genes <- unique(paad_deps$hugo_symbol)
jun_targets <- jun_bs$gene


fisher_tf_in_paad <- function(tf_genes, alternative = "two.sided") {
  tf_genes <- unlist(unique(tf_genes))
  tf_genes <- tf_genes[tf_genes %in% all_paad_genes]

  if (length(tf_genes) < 5) {
    return(NULL)
  }

  A <- length(setdiff(all_paad_genes, c(tf_genes, sig_paad_genes)))
  B <- length(setdiff(sig_paad_genes, tf_genes))
  C <- length(setdiff(tf_genes, sig_paad_genes))
  D <- length(intersect(sig_paad_genes, tf_genes))

  stopifnot(sum(c(A, B, C, D)) == n_distinct(all_paad_genes))

  m <- matrix(c(A, B, C, D), nrow = 2, byrow = FALSE)
  fish_res <- fisher.test(m, alternative = alternative)
  return(fish_res)
}


tf_stats <- encode_tf_bindingsites %>%
  group_by(gene_set) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    fish_res = purrr::map(data, fisher_tf_in_paad),
    fish_tidy = purrr::map(fish_res, ~ tidy(.x))
  ) %>%
  unnest(fish_tidy) %>%
  janitor::clean_names()

tf_stats %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
  filter(p_value < 0.05 & p_value_adj < 0.25)

if (nrow(tf_stats) > 0) {
  tf_stats_volcano <- tf_stats %>%
    mutate(
      tf_gs_len = map_dbl(data, ~ length(unlist(.x))),
      norm_log_or = log2(estimate) / sqrt(tf_gs_len),
      jun = ifelse(gene_set == "JUN", "JUN", "other"),
      labels = ifelse(p_value < 0.05 | gene_set == "JUN", gene_set, NA)
    ) %>%
    ggplot(aes(x = log2(estimate), y = -log10(p_value))) +
    geom_hline(yintercept = 0, size = 0.5, color = "grey25") +
    geom_vline(xintercept = 0, size = 0.5, color = "grey25") +
    geom_hline(
      yintercept = -log10(0.05),
      size = 0.7, linetype = 2, color = "tomato"
    ) +
    geom_point(aes(color = jun, size = jun), alpha = 0.6) +
    ggrepel::geom_text_repel(aes(label = labels),
      size = 2, color = "black",
      family = "Arial"
    ) +
    scale_color_manual(
      values = c("dodgerblue", "grey50"),
      guide = FALSE
    ) +
    scale_size_manual(values = c(0.9, 0.7), guide = FALSE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "log2 OR", y = "-log10 p-value")
  ggsave_wrapper(
    tf_stats_volcano,
    plot_path(GRAPHS_DIR, "tf_stats_volcano.svg"),
    "medium"
  )
}


#### ---- JUN vs MAPK8 gene effect ---- ####

# Save scatter plot for Fig.
fignum <- 14
supp <- TRUE

# JUN vs. MAPK8 gene effect

JUN_MAPK8_df <- junpw_depmap %>%
  filter(hugo_symbol %in% c("JUN", "MAPK8")) %>%
  select(
    dep_map_id, kras_allele, hugo_symbol,
    gene_effect, rna_expression
  ) %>%
  pivot_wider(c(dep_map_id, kras_allele),
    names_from = hugo_symbol,
    values_from = c(rna_expression, gene_effect)
  )

jun_mapk8_fit <- lm(gene_effect_MAPK8 ~ gene_effect_JUN, data = JUN_MAPK8_df)
summary(jun_mapk8_fit)

fit_x <- round(coef(jun_mapk8_fit)[[2]], 3)
fit_b <- round(coef(jun_mapk8_fit)[[1]], 3)

jun_mapk8_model_df <- tibble(
  x = -0.62,
  y = c(0.52, 0.49, 0.46),
  label = c(
    glue("*y* = {fit_x}*x* + {fit_b}"),
    glue("p-value: {round(glance(jun_mapk8_fit)$p.value, 3)}"),
    glue("*R*<sup>2</sup>: {round(glance(jun_mapk8_fit)$r.squared, 2)}")
  )
)

JUN_MAPK8_scatter <- JUN_MAPK8_df %>%
  ggplot(aes(x = gene_effect_JUN, gene_effect_MAPK8)) +
  geom_smooth(method = "lm", color = "grey25", size = 0.8, linetype = 2) +
  geom_point(aes(color = kras_allele)) +
  geom_richtext(aes(x = x, y = y, label = label),
    data = jun_mapk8_model_df,
    size = 2.3, family = "arial", hjust = 0,
    fill = NA, label.color = NA,
    label.padding = grid::unit(rep(0, 4), "pt")
  ) +
  scale_color_manual(values = short_allele_pal) +
  scale_shape_manual(values = c(17, 16)) +
  scale_size_manual(values = c(1, 2, 3)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank())
ggsave_wrapper(
  JUN_MAPK8_scatter,
  plot_path(GRAPHS_DIR, "JUN_MAPK8_scatter.svg"),
  "small"
)
saveFigRds(JUN_MAPK8_scatter, "JUN_MAPK8_scatter")
