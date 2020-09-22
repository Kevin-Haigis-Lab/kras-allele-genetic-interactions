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


#### ---- Comparison with allele-specific comutation analysis ---- ####


get_allelespecific_genes <- function(cancer, genetic_interaction) {
  genetic_interaction_df %>%
    filter(
      cancer == !!cancer & genetic_interaction == !!genetic_interaction
    ) %>%
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
    count_size = 3
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
    labs(fill = "comutation analysis",
         title = title)
}


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
    plot_path(GRAPHS_DIR, glue("comutation-venn-patch_{cancer}.svg")),
    "wide"
  )
})




