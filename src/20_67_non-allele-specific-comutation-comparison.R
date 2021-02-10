# Analyze the results of the non-allele-specific comutation analysis.

library(mustashe)

GRAPHS_DIR <- "20_67_non-allele-specific-comutation-comparison"
reset_graph_directory(GRAPHS_DIR)

source(file.path("src", "20_34_rc-fisher-assessment-processes.R"))

theme_set(theme_bw(base_size = 7, base_family = "Arial"))


#### ---- Prepare non-allele-specific RC test results ---- ####

nonallele_rc_res <- rc_test_nonallele_results %>%
  select(-rc_test_type, -allele) %>%
  rename(rc_test_type = kras_allele) %>%
  select(
    cancer, hugo_symbol, p_val, t_AM, t_BM_ge,
    num_samples_per_cancer:num_mut_per_cancer_allele
  ) %>%
  assess_rc_test_significance()


#### ---- Prepare non-allele-specific Fisher test results ---- ####

#             KRAS
#
#             WT  M
#           +---+---+
#        WT | A | C |
#   gene    +-------+
#         M | B | D |
#           +---+---+
#
# Parse a contingency table of the above format into `n00` -> `n11` columns.
parse_comutation_contingency_table <- function(tbl) {
  tibble(
    n00 = unlist(tbl[1, 1]),
    n01 = unlist(tbl[1, 2]),
    n10 = unlist(tbl[2, 1]),
    n11 = unlist(tbl[2, 2])
  )
}

stash(
  "nonallele_fish_res",
  depends_on = "nonallele_specific_increased_comutation_df",
  {
    nonallele_fish_res <- nonallele_specific_increased_comutation_df %>%
      filter(hugo_symbol != "KRAS") %>%
      mutate(
        comut_tib = map(comut_ct_tbl, parse_comutation_contingency_table),
        p_value_less = 1
      ) %>%
      rename(p_value_great = p_value) %>%
      unnest(comut_tib) %>%
      assess_fisher_test_significance()
  }
)


#### ---- Venn diagrams ---- ####

get_allelespecific_genes <- function(cancer,
                                     genetic_interaction,
                                     allele = "all") {
  genetic_interaction_df %>%
    filter(
      cancer == !!cancer & genetic_interaction == !!genetic_interaction
    ) %>%
    filter(allele %in% !!allele | !!allele == "all") %>%
    u_pull(hugo_symbol)
}


get_nonallelespecific_genes <- function(cancer, genetic_interaction) {
  if (genetic_interaction == "comutation") {
    gs <- nonallele_fish_res %>%
      filter(is_sig & cancer == !!cancer) %>%
      u_pull(hugo_symbol)
  } else {
    gs <- nonallele_rc_res %>%
      filter(cancer == !!cancer & is_sig) %>%
      u_pull(hugo_symbol)
  }

  return(gs)
}


make_comutation_venn <- function(cancer, genetic_interaction, title = NULL) {
  fill_pal <- c(
    "allele-specific" = "orange",
    "non-allele-specific" = "blue"
  )

  ggvenndiagram(
    s1 = get_allelespecific_genes(cancer, genetic_interaction),
    s2 = get_nonallelespecific_genes(cancer, genetic_interaction),
    circle_fill = c("allele-specific", "non-allele-specific"),
    circle_alpha = 0.2,
    cat_family = "Arial",
    count_size = 2.7
  ) +
    scale_fill_manual(
      values = fill_pal,
      guide = guide_legend(override.aes = list(size = 0, alpha = 0.6))
    ) +
    coord_fixed() +
    theme_void(base_size = 7, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.key.size = unit(3, "mm"),
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 7)
    ) +
    labs(
      fill = "comutation analysis",
      title = title
    )
}

# Save for Supp. Data.
nonallele_rc_res %>%
  filter(is_sig) %>%
  select(-is_sig) %>%
  save_supp_data(
    num = 11,
    sheet_name = "Non-allele-specific comutation"
  )


walk(unique(nonallele_rc_res$cancer), function(cancer) {
  rc_venn <- make_comutation_venn(
    cancer,
    "exclusivity",
    "Reduced comutation"
  )
  fish_venn <- make_comutation_venn(
    cancer,
    "comutation",
    "Increased comutation"
  )

  patch_venn <- ((rc_venn | fish_venn) / guide_area()) +
    plot_layout(
      guides = "collect",
      heights = c(4, 1)
    ) +
    plot_annotation(
      title = cancer,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )

  ggsave_wrapper(
    patch_venn,
    plot_path(
      GRAPHS_DIR,
      as.character(glue("comutation-venn-patch_{cancer}.svg"))
    ),
    "wide"
  )

  saveFigRds(rc_venn, as.character(glue("rc-venn-diagram_{cancer}")))
  saveFigRds(fish_venn, as.character(glue("fish-venn-diagram_{cancer}")))
})




#### ---- Bar-plots for number of new and lost genes per KRAS allele ---- ####



genetic_interaction_sets <- genetic_interaction_df %>%
  distinct(cancer, allele, genetic_interaction)

genetic_interaction_sets$allele_sets <- pmap(
  genetic_interaction_sets,
  get_allelespecific_genes
)

genetic_interaction_sets %<>%
  mutate(
    cancer_sets = map2(
      cancer,
      genetic_interaction,
      get_nonallelespecific_genes
    ),
    new_genes = map2_dbl(
      allele_sets, cancer_sets, ~ length(setdiff(.x, .y))
    ),
    old_genes = map2_dbl(
      allele_sets, cancer_sets, ~ length(intersect(.x, .y))
    )
  )


comutation_comparison_barplot <- function(cancer,
                                          data,
                                          y,
                                          title = NULL,
                                          fn_glue = NA) {
  pos <- position_dodge(width = 0.9)
  p <- data %>%
    ggplot(aes(x = allele, y = {{ y }}, fill = comutation)) +
    facet_wrap(~which_genes_label, scales = "free_y", nrow = 1) +
    geom_col(position = pos) +
    geom_text(
      aes(y = is_tested_label_y, label = is_tested_label),
      position = pos,
      color = "grey20",
      size = 1.8,
      family = "Arial"
    ) +
    scale_fill_manual(values = comut_updown_pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      legend.key.size = unit(3, "mm")
    ) +
    labs(
      x = "KRAS allele",
      y = "number of genes",
      fill = "comutation\ninteraction",
      title = title
    )

  if (!is.na(fn_glue)) {
    fn <- as.character(glue(fn_glue))
    ggsave_wrapper(
      p,
      plot_path(GRAPHS_DIR, fn),
      "wide"
    )
    saveFigRds(p, fn)
  }
  invisible(p)
}


gene_comparison_names <- tribble(
  ~which_genes, ~which_genes_label,
  "new_genes", "Only found in the allele-specific analysis",
  "old_genes", "Found in the non-allele-specific analysis"
)

genetic_interaction_sets %>%
  mutate(
    allele = factor_alleles(allele),
    comutation = switch_comut_terms(genetic_interaction),
    is_tested_label = ""
  ) %>%
  filter(!(cancer == "MM" & allele == "Q61L")) %>%
  complete(
    nesting(cancer, allele),
    comutation,
    fill = list(new_genes = 0, old_genes = 0, is_tested_label = "NA")
  ) %>%
  select(-allele_sets, -cancer_sets) %>%
  pivot_longer(
    c(new_genes, old_genes),
    names_to = "which_genes",
    values_to = "n_genes"
  ) %>%
  group_by(cancer, which_genes) %>%
  mutate(is_tested_label_y = n_genes + (0.04 * max(n_genes))) %>%
  ungroup() %>%
  left_join(gene_comparison_names, by = "which_genes") %>%
  group_by(cancer) %>%
  nest() %>%
  ungroup() %>%
  pwalk(
    comutation_comparison_barplot,
    y = n_genes,
    fn_glue = "comutation-comparison_bar_{cancer}.svg"
  )
