# Analysis of the predicted KRAS alleles per sample.

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"
reset_graph_directory(GRAPHS_DIR)

theme_set(theme_bw(base_size = 7, base_family = "Arial"))

# Filter out infrequent KRAS alleles
filter_kras_allele_ct <- function(df, min_n = 20) {
    original_groups <- dplyr::groups(df)
    df %>%
        group_by(cancer, real_kras_allele) %>%
        filter(n_distinct(tumor_sample_barcode) >= !!min_n) %>%
        group_by(!!!original_groups)
}


real_kras_mutations <- trinucleotide_mutations_df %>%
  distinct(cancer, tumor_sample_barcode, kras_allele)

ranked_allele_predictions <- kras_allele_predictions %>%
    inner_join(real_kras_mutations %>% rename(real_kras_allele = kras_allele),
        by = c("cancer", "tumor_sample_barcode")
    ) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    arrange(-allele_prob) %>%
    mutate(allele_idx = row_number()) %>%
    ungroup() %>%
    arrange(cancer, tumor_sample_barcode, allele_idx)



allele_prob_ranking_barplot <- ranked_allele_predictions %>%
  filter(real_kras_allele == kras_allele) %>%
  count(cancer, real_kras_allele, allele_idx) %>%
  mutate(
    top_2 = allele_idx < 3,
    correct_allele = allele_idx == 1,
    real_kras_allele = fct_reorder(real_kras_allele,
      -correct_allele,
      .fun = mean
    )
  ) %>%
  ggplot(aes(real_kras_allele, n)) +
  facet_wrap(~cancer, nrow = 2, scales = "free") +
  geom_col(aes(fill = fct_rev(factor(allele_idx)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(type = "div", palette = "RdYlBu") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    x = NULL,
    y = "number of tumor samples",
    fill = "rank",
    subtitle = "For each cancer, the alleles are ordered by the accuracy of the mutational signature prediction",
    title = "The rank of the observed allele when ordered by probability using mutational signatures"
  )
ggsave_wrapper(
  allele_prob_ranking_barplot,
  plot_path(GRAPHS_DIR, "allele-prob-ranking-barplot.svg"),
  "medium"
)



ranked_probability_plot <- function(df, allele, cancer, ignore_alleles = c()) {
  title <- glue("The distribution of probabilities of each KRAS allele in {cancer} tumors with KRAS {allele}")
  df %>%
    filter(real_kras_allele == !!allele & cancer %in% !!cancer) %>%
    filter(!(kras_allele %in% !!ignore_alleles)) %>%
    mutate(kras_allele = fct_reorder(kras_allele, allele_prob)) %>%
    group_by(tumor_sample_barcode) %>%
    mutate(hit = any(allele_idx == 1 & kras_allele == real_kras_allele)) %>%
    ungroup() %>%
    mutate(hit = factor(hit, levels = c("TRUE", "FALSE"))) %>%
    ggplot(aes(x = kras_allele, y = allele_prob)) +
    geom_line(aes(group = tumor_sample_barcode, color = hit), alpha = 0.3) +
    geom_point(aes(color = hit), alpha = 1, size = 2) +
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.02)),
      limits = c(0, NA)
    ) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(
      x = "possible KRAS allele",
      y = "probability of allele from mutational signatures",
      color = glue("correct prediction"),
      title = title
    )
}


save_ranked_probability_plot <- function(cancer, allele, plt, ...) {
  fname <- glue("ranked-probabilit-lines_{cancer}_{allele}.svg")
  ggsave_wrapper(
    plt,
    plot_path(GRAPHS_DIR, fname),
    "wide"
  )
}

ranked_prob_plot_data <- ranked_allele_predictions %>%
    filter_kras_allele_ct() %>%
    distinct(cancer, real_kras_allele) %>%
    left_join(
        tibble(
            real_kras_allele = c("G12C"),
            ignore_alleles = c("G13C")
        ),
        by = "real_kras_allele"
    ) %>%
    rename(allele = real_kras_allele)

ranked_prob_plot_data$plt <- pmap(
    ranked_prob_plot_data,
    ranked_probability_plot,
    df = ranked_allele_predictions
)

pwalk(ranked_prob_plot_data, save_ranked_probability_plot)




allele_predictions_graph_plot <- function(gr, cancer) {
  gr_arrow <- arrow(angle = 25, length = unit(1.6, "mm"), type = "closed")
  end_cap <- circle(4, "mm")
  gr %E>%
    filter(cancer == !!cancer) %N>%
    jhcutils::filter_component_size(min_size = 2) %>%
    ggraph(layout = "circle") +
    geom_edge_parallel(aes(width = total_prob, color = edge_color),
      alpha = 0.9,
      arrow = gr_arrow,
      end_cap = end_cap
    ) +
    geom_edge_loop(aes(width = total_prob, color = edge_color),
      alpha = 0.9,
      arrow = gr_arrow,
      end_cap = end_cap
    ) +
    geom_node_label(aes(label = name, fill = name, color = lbl_color),
      label.size = 0, fontface = "bold",
      size = 2,
      family = "Arial"
    ) +
    scale_edge_width_continuous(range = c(0.3, 1.6)) +
    scale_edge_color_manual(values = short_allele_pal, drop = TRUE, guide = FALSE) +
    scale_fill_manual(values = short_allele_pal, drop = TRUE, guide = FALSE) +
    scale_color_manual(values = c("black", "white"), guide = FALSE) +
    theme(
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    labs(
      edge_width = "mutational signature\nestimated probability",
      title = glue("The average probability of {cancer} tumor samples\nwith one KRAS allele to be predicted to have another")
    )
}


ranked_allele_predictions_gr <- ranked_allele_predictions %>%
  group_by(cancer) %>%
  filter(real_kras_allele %in% c("WT", kras_allele)) %>%
  group_by(cancer, kras_allele, real_kras_allele) %>%
  summarise(total_prob = mean(allele_prob)) %>%
  ungroup() %>%
  rename(
    from = real_kras_allele,
    to = kras_allele
  ) %>%
  as_tbl_graph(directed = TRUE) %E>%
  mutate(edge_color = .N()$name[from]) %N>%
  mutate(lbl_color = name %in% kras_dark_lbls) %>%
  jhcutils::filter_component_size(min_size = 2)

for (cancer in unique(ranked_allele_predictions$cancer)) {
  fn <- glue("predicted-allele-graph_{cancer}.svg")
  p <- allele_predictions_graph_plot(
    gr = ranked_allele_predictions_gr,
    cancer = cancer
  )
  ggsave_wrapper(p, plot_path(GRAPHS_DIR, fn), "small")
}



most_likely_other_allele <- ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  filter(allele_idx == 1) %>%
  count(cancer, real_kras_allele, kras_allele) %>%
  mutate(
    real_kras_allele = factor_alleles(real_kras_allele),
    kras_allele = factor_alleles(kras_allele)
  ) %>%
  ggplot(aes(x = real_kras_allele, y = n)) +
  facet_wrap(~cancer, scales = "free", nrow = 2) +
  geom_col(aes(fill = kras_allele), position = "fill", width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = short_allele_pal) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "observed KRAS allele",
    y = "fraction of tumor samples most most likely to be another allele",
    fill = "most likely\nKRAS allele",
    legend.key.size = unit(3, "mm")
  )
ggsave_wrapper(
  most_likely_other_allele,
  plot_path(GRAPHS_DIR, "most-likely-other-allele.svg"),
  "wide"
)


prob_of_other_alleles <- ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  group_by(cancer, real_kras_allele, kras_allele) %>%
  summarise(total_prob = sum(allele_prob)) %>%
  ungroup() %>%
  mutate(
    real_kras_allele = factor_alleles(real_kras_allele),
    kras_allele = factor_alleles(kras_allele)
  ) %>%
  ggplot(aes(x = real_kras_allele, y = total_prob)) +
  facet_wrap(~cancer, scales = "free", nrow = 2) +
  geom_col(aes(fill = kras_allele), position = "fill", width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = short_allele_pal) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "observed KRAS allele",
    y = "fraction of tumor samples most most likely to be another allele",
    fill = "alternative\nKRAS allele",
    legend.key.size = unit(3, "mm")
  )
ggsave_wrapper(
  prob_of_other_alleles,
  plot_path(GRAPHS_DIR, "cum-prob-of-other-allele.svg"),
  "wide"
)


distribution_of_allele_probs <- ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  mutate(correct_allele = kras_allele == real_kras_allele) %>%
  ggplot(aes(x = allele_prob)) +
  facet_grid(real_kras_allele ~ cancer, scales = "free") +
  geom_density(aes(
    color = kras_allele, fill = kras_allele,
    alpha = correct_allele
  )) +
  scale_color_manual(values = short_allele_pal, drop = TRUE) +
  scale_fill_manual(values = short_allele_pal, drop = TRUE) +
  scale_alpha_manual(values = c(0, 0.25), guide = FALSE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0.25, 1, 0.25)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    x = "probability of KRAS allele from mutational signatures",
    y = "probability density",
    color = "KRAS allele",
    fill = "KRAS allele",
    title = "Distributions of probabilities of KRAS alleles",
    subtitle = "The tumor samples are separated vertically by the KRAS mutation they contain."
  )
ggsave_wrapper(
  distribution_of_allele_probs,
  plot_path(GRAPHS_DIR, "allele-prob-distribution.svg"),
  "large"
)
