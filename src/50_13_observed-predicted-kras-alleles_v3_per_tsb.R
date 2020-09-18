# Analysis of the predicted KRAS alleles per sample.

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"
reset_graph_directory(GRAPHS_DIR)

theme_set(theme_bw(base_size = 7, base_family = "Arial"))

# Minumum number of tumor samples to include the KRAS allele in the analysis.
MIN_ALLELE_TSB <- 20

filter_kras_allele_ct <- function(df) {
  df %>%
    filter(num_allele_tsb >= MIN_ALLELE_TSB)
}


real_kras_mutations <- trinucleotide_mutations_df %>%
  distinct(cancer, tumor_sample_barcode, kras_allele)


prepare_ranked_allele_prediction_df <- function(df) {
  df %>%
    inner_join(
      real_kras_mutations %>% rename(real_kras_allele = kras_allele),
      by = c("tumor_sample_barcode", "cancer")
    ) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    arrange(-allele_prob) %>%
    mutate(allele_idx = row_number()) %>%
    group_by(cancer, real_kras_allele) %>%
    mutate(num_allele_tsb = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    arrange(cancer, tumor_sample_barcode, allele_idx)
}


ranked_allele_predictions <- kras_allele_predictions %>%
  prepare_ranked_allele_prediction_df()

all_ranked_allele_predictions <- all_kras_allele_predictions %>%
  prepare_ranked_allele_prediction_df()



#### ---- Bar-plots of ranking of predictions ---- ####


plot_allele_probability_barplot <- function(df,
                                            geombar_position = "stack",
                                            max_pred_idx = 11) {

  if (max(df$allele_idx) > max_pred_idx) {
    fill_idx_lbls <- rev(c(
      as.character(seq(1, max_pred_idx - 1)),
      glue("≥{max_pred_idx}")
    ))
  } else {
    fill_idx_lbls <- as.character(seq(1, max(df$allele_idx)))
  }

  df %>%
    filter(real_kras_allele == kras_allele) %>%
    count(cancer, real_kras_allele, allele_idx) %>%
    mutate(
      top_2 = allele_idx < 3,
      correct_allele = allele_idx == 1,
      real_kras_allele = fct_reorder(real_kras_allele,
        -correct_allele,
        .fun = mean
      ),
      allele_idx = minmax(allele_idx, 0, 11)
    ) %>%
    ggplot(aes(real_kras_allele, n)) +
    facet_wrap(~cancer, nrow = 2, scales = "free") +
    geom_col(
      aes(fill = fct_rev(factor(allele_idx))),
      position = geombar_position
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_fill_brewer(
      type = "div",
      palette = "RdYlBu",
      labels = fill_idx_lbls
    ) +
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
}

ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  plot_allele_probability_barplot() %>%
  ggsave_wrapper(
    plot_path(GRAPHS_DIR, "allele-prob-ranking_stack_barplot.svg"),
    "medium"
  )

p <- ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  plot_allele_probability_barplot(geombar_position = "fill")
p <- p +
  scale_y_continuous(expand = c(0, 0))
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "allele-prob-ranking_fill_barplot.svg"),
  "medium"
)


all_ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  plot_allele_probability_barplot(max_pred_idx = 7) %>%
  ggsave_wrapper(
    plot_path(GRAPHS_DIR, "all-allele-prob-ranking_stack_barplot.svg"),
    "medium"
  )

p <- all_ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  plot_allele_probability_barplot(geombar_position = "fill", max_pred_idx = 7)
p <- p +
  scale_y_continuous(expand = c(0, 0))
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "all-allele-prob-ranking_fill_barplot.svg"),
  "medium"
)



#### ---- Line-plots for tumor samples of a single KRAS allele ---- ####


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



#### ---- Graphs of prob. of other alleles for each cancer type ---- ####

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



#### ---- Bar-plots of probability of other alleles ---- ####


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



#### ---- Distributions of allele probabilities ---- ####

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



allele_prob_per_allele_bars <- ranked_allele_predictions %>%
  filter_kras_allele_ct(min_n = 10) %>%
  group_by(cancer) %>%
  filter(kras_allele %in% real_kras_allele) %>%
  group_by(cancer, real_kras_allele) %>%
  mutate(num_tsb = n_distinct(tumor_sample_barcode)) %>%
  group_by(cancer, real_kras_allele, kras_allele) %>%
  summarise(total_prob = sum(allele_prob) / unique(num_tsb)) %>%
  ungroup() %>%
  mutate(
    correct_allele = real_kras_allele == kras_allele,
    real_kras_allele = factor_alleles(real_kras_allele),
    kras_allele = factor_alleles(kras_allele)
  ) %>%
  ggplot(aes(x = real_kras_allele, y = total_prob)) +
  facet_wrap(~cancer, scales = "free", nrow = 2) +
  geom_col(
    aes(fill = kras_allele, alpha = correct_allele),
    position = "dodge", width = 0.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = short_allele_pal) +
  scale_alpha_manual(values = c(0.5, 1), guide = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "observed KRAS allele",
    y = "average probability of KRAS allele",
    fill = "alternative\nKRAS allele",
    legend.key.size = unit(3, "mm")
  )
ggsave_wrapper(
  allele_prob_per_allele_bars,
  plot_path(GRAPHS_DIR, "allele-prob-per-allele-bars.svg"),
  "medium"
)


allele_prob_distribution_per_kras_mutation <- function(df, cancer, allele) {
  xlbl <- glue("probability of {allele} allele")
  title <- glue("Probability of a {allele} allele in {cancer} tumors\nwith various KRAS alleles")

  df %>%
    filter(cancer == !!cancer & kras_allele == !!allele) %>%
    mutate(real_kras_allele = fct_lump_min(real_kras_allele, min = 20)) %>%
    ggplot(aes(allele_prob)) +
    geom_density(
      aes(color = real_kras_allele, fill = real_kras_allele),
      alpha = 0.2,
      adjust = 1.5
    ) +
    scale_color_manual(values = short_allele_pal, drop = TRUE) +
    scale_fill_manual(values = short_allele_pal, drop = TRUE) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme(
      axis.ticks = element_blank(),
      legend.position = c(0.7, 0.8),
      legend.background = element_rect(fill = "white", color = "grey50"),
      legend.key.size = unit(3, "mm")
    ) +
    labs(
      x = xlbl,
      y = "density",
      color = "observed KRAS allele",
      fill = "observed KRAS allele",
      title = title
    )
}


allele_prob_dist_data <- ranked_allele_predictions %>%
  distinct(cancer, kras_allele) %>%
  arrange(cancer, kras_allele) %>%
  rename(allele = kras_allele)

allele_prob_dist_data$plt <- pmap(
  allele_prob_dist_data,
  allele_prob_distribution_per_kras_mutation,
  df = ranked_allele_predictions
)

pwalk(allele_prob_dist_data, function(cancer, allele, plt, ...) {
  ggsave_wrapper(
    plt,
    plot_path(
      GRAPHS_DIR,
      glue("allele-prob-other-kras-tumors_{cancer}_{allele}.svg")
    ),
    "small"
  )
})




average_allele_prob_incorrect_alleles <- ranked_allele_predictions %>%
  filter(kras_allele != real_kras_allele) %>%
  group_by(cancer, kras_allele) %>%
  summarise(avg_incorrect_allele_prob = mean(allele_prob)) %>%
  ungroup()

diff_between_obs_and_other_plt <- ranked_allele_predictions %>%
  filter_kras_allele_ct() %>%
  filter(kras_allele == real_kras_allele) %>%
  select(-kras_allele) %>%
  left_join(average_allele_prob_incorrect_alleles,
    by = c("cancer", "real_kras_allele" = "kras_allele")
  ) %>%
  mutate(adj_allele_prob = allele_prob - avg_incorrect_allele_prob) %>%
  ggplot(aes(x = adj_allele_prob)) +
  facet_wrap(~cancer, scales = "free", nrow = 2) +
  geom_vline(xintercept = 0, lty = 2, size = 0.6, color = "grey20") +
  geom_density(
    aes(color = real_kras_allele, fill = real_kras_allele),
    alpha = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = short_allele_pal, drop = TRUE) +
  scale_fill_manual(values = short_allele_pal, drop = TRUE) +
  theme(
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "probability of observed KRAS allele - average probability in tumors with other KRAS allele",
    y = "density",
    title = "The distribution of the difference in probability of a KRAS allele between\nsamples with the KRAS allele and those with another allele"
  )

ggsave_wrapper(
  diff_between_obs_and_other_plt,
  plot_path(GRAPHS_DIR, "diff-between-obs-and-other-density.svg"),
  "medium"
)


# Find the fraction of `non_mut_df$allele_prob` less than `prob`.
calculate_frac_with_lower_prob <- function(prob, non_mut_df) {
  mean(non_mut_df$allele_prob < prob)
}


# Plot the distribution of allele probabilities for one allele in tumors with
# the allele vs. those without.
allele_prob_dist_plot1 <- function(df, allele) {
  pal <- c(
    "TRUE" = short_allele_pal[[allele]],
    "FALSE" = "grey50"
  )
  lbls <- c(
    "TRUE" = glue("{allele} mutant"),
    "FALSE" = glue("non-{allele} mutant")
  )
  xlbl <- glue("probability of {allele} allele")

  df %>%
    ggplot(aes(x = allele_prob)) +
    geom_density(
      aes(color = is_allele, fill = is_allele),
      alpha = 0.5,
      adjust = 1.3
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_color_manual(values = pal, labels = lbls) +
    scale_fill_manual(values = pal, labels = lbls) +
    theme(
      axis.ticks = element_blank(),
      legend.position = c(0.75, 0.8),
      legend.key.size = unit(3, "mm")
    ) +
    labs(
      x = xlbl,
      y = "density",
      color = NULL,
      fill = NULL
    )
}


# Plot the distribution of the fraction of probabilities with lower allele_prob
# without the allele.
allele_prob_frac_dist_plot2 <- function(df, allele) {
  xlbl <- glue("fraction of non-{allele} mutant samples\nwith a lower probability of {allele} mutation")
  df %>%
    ggplot(aes(x = frac_less_prob)) +
    geom_vline(xintercept = 0.5, lty = 2, color = "grey20", size = 0.6) +
    geom_density(fill = "#8491a1", alpha = 0.5, size = 1, color = "#8491a1") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      x = xlbl,
      y = "density"
    )
}


get_allele_prob_lesser_frac <- function(cancer, allele, prob_data,
                                        ignore_alleles = NULL) {
  prob_allele <- prob_data %>%
    filter(cancer == !!cancer) %>%
    filter(kras_allele == !!allele & real_kras_allele == !!allele)

  prob_not_allele <- prob_data %>%
    filter(cancer == !!cancer) %>%
    filter(kras_allele == !!allele & real_kras_allele != !!allele) %>%
    filter(!(kras_allele %in% !!ignore_alleles))

  prob_allele %>%
    mutate(frac_less_prob = map_dbl(allele_prob,
      calculate_frac_with_lower_prob,
      non_mut_df = prob_not_allele
    ))
}


fraction_of_allele_prob_plot <- function(cancer,
                                         allele,
                                         prob_data,
                                         ignore_alleles = NULL) {
  p1 <- prob_data %>%
    filter(cancer == !!cancer & kras_allele == !!allele) %>%
    mutate(is_allele = real_kras_allele == !!allele) %>%
    allele_prob_dist_plot1(allele = allele)


  p2 <- get_allele_prob_lesser_frac(
    cancer = cancer,
    allele = allele,
    prob_data = prob_data,
    ignore_alleles = ignore_alleles
  ) %>%
    allele_prob_frac_dist_plot2(allele = allele)

  patch_title <- glue("The probability of {allele} mutations in {cancer} tumors with and without a {allele} allele")
  patch <- (p1 | p2) +
    plot_annotation(
      title = patch_title,
      theme = theme(plot.title = element_text(hjust = 0.5))
    )

  fn <- glue("allele-prob-fraction-greater_{cancer}_{allele}.svg")
  ggsave_wrapper(
    patch,
    plot_path(GRAPHS_DIR, fn),
    "wide"
  )
}


ranked_allele_predictions %>%
  group_by(cancer) %>%
  filter(real_kras_allele %in% kras_allele) %>%
  ungroup() %>%
  distinct(cancer, real_kras_allele) %>%
  rename(allele = real_kras_allele) %>%
  left_join(tibble(
    cancer = "LUAD",
    allele = c("G12C", "G13C"),
    ignore_alleles = c("G13C", "G12C")
  ),
  by = c("cancer", "allele")
  ) %>%
  pwalk(fraction_of_allele_prob_plot, prob_data = ranked_allele_predictions)





all_dist_frac_lower_prob <- ranked_allele_predictions %>%
  group_by(cancer) %>%
  filter(real_kras_allele %in% kras_allele) %>%
  ungroup() %>%
  distinct(cancer, real_kras_allele) %>%
  rename(allele = real_kras_allele) %>%
  left_join(tibble(
    cancer = "LUAD",
    allele = c("G12C", "G13C"),
    ignore_alleles = c("G13C", "G12C")
  ),
  by = c("cancer", "allele")
  ) %>%
  pmap(get_allele_prob_lesser_frac,
    prob_data = ranked_allele_predictions
  ) %>%
  bind_rows() %>%
  ggplot(aes(frac_less_prob)) +
  facet_grid(kras_allele ~ cancer, scales = "free_y") +
  geom_vline(
    xintercept = 0.5,
    lty = 2,
    size = 0.6,
    color = "grey25"
  ) +
  geom_density(
    aes(color = kras_allele, fill = kras_allele),
    alpha = 0.6
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(
    values = short_allele_pal,
    drop = TRUE
  ) +
  scale_fill_manual(
    values = short_allele_pal,
    drop = TRUE
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    x = "fraction of tumors samples without the allele mutant with a lower probability of the allele",
    y = "density",
    color = "KRAS allele",
    fill = "KRAS allele"
  )

ggsave_wrapper(
  all_dist_frac_lower_prob,
  plot_path(GRAPHS_DIR, "all-dist-frac-lower-prob.svg"),
  "medium"
)



#### ---- All KRAS alleles ---- ####





max_pred_idx <- 7
fill_idx_lbls <- rev(c(
  as.character(seq(1, max_pred_idx - 1)),
  glue("≥{max_pred_idx}")
))

all_ranked_allele_predictions %>%
  filter(num_allele_tsb > 10) %>%
  filter(real_kras_allele == kras_allele) %>%
  count(cancer, real_kras_allele, allele_idx) %>%
  mutate(
    top_2 = allele_idx < 3,
    correct_allele = allele_idx == 1,
    real_kras_allele = fct_reorder(real_kras_allele,
      -correct_allele,
      .fun = mean
    ),
    allele_idx = minmax(allele_idx, 0, 7)
  ) %>%
  ggplot(aes(real_kras_allele, n)) +
  facet_wrap(~cancer, nrow = 2, scales = "free") +
  geom_col(aes(fill = fct_rev(factor(allele_idx))), position = "stack") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(type = "div", palette = "RdYlBu", labels = fill_idx_lbls) +
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



prob_all_alleles <- all_ranked_allele_predictions %>%
  group_by(cancer, real_kras_allele, kras_allele, num_allele_tsb) %>%
  summarise(cum_allele_prob = sum(allele_prob) / unique(num_allele_tsb)) %>%
  ungroup() %>%
  filter(num_allele_tsb > 10) %>%
  ggplot(aes(real_kras_allele, cum_allele_prob)) +
  facet_wrap(~cancer, nrow = 2, scales = "free_x") +
  geom_col(aes(fill = kras_allele)) +
  scale_fill_manual(
    values = short_allele_pal,
    drop = TRUE
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.ticks = element_blank(),
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    x = "observed KRAS allele",
    y = "average probability of the alternative KRAS allele",
    fill = "alternative\nKRAS allele"
  )
ggsave_wrapper(
  prob_all_alleles,
  plot_path(GRAPHS_DIR, "prob-all-alleles.svg"),
  "medium"
)



all_ranked_allele_predictions %>%
  filter(cancer == "COAD" & real_kras_allele == "G12D") %>%
  mutate(
    correct_allele = real_kras_allele == kras_allele,
    kras_allele = fct_reorder(kras_allele, allele_prob)
  ) %>%
  ggplot(aes(x = allele_prob)) +
  geom_density_ridges(
    aes(
      y = kras_allele,
      color = kras_allele,
      fill = kras_allele,
      alpha = correct_allele
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_color_manual(
    values = short_allele_pal,
    drop = TRUE,
    guide = FALSE
  ) +
  scale_fill_manual(
    values = short_allele_pal,
    drop = TRUE,
    guide = FALSE
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 0.6, "FALSE" = 0.2),
    guide = FALSE
  ) +
  theme(
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    x = "probability of allele from mutational signatures",
    y = "possible KRAS allele",
    title = "Probabilities of common KRAS mutations in COAD tumors samples with G12D"
  )


all_ranked_allele_predictions %>%
  filter(cancer == "COAD") %>%
  filter(kras_allele == "G12D") %>%
  mutate(
    correct_allele = real_kras_allele == kras_allele,
    real_kras_allele = fct_reorder(real_kras_allele, allele_prob)
  ) %>%
  filter(num_allele_tsb > 10) %>%
  ggplot(aes(x = allele_prob)) +
  geom_density_ridges(
    aes(
      y = real_kras_allele,
      color = real_kras_allele,
      fill = real_kras_allele,
      alpha = correct_allele
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_color_manual(
    values = short_allele_pal,
    drop = TRUE,
    guide = FALSE
  ) +
  scale_fill_manual(
    values = short_allele_pal,
    drop = TRUE,
    guide = FALSE
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 0.6, "FALSE" = 0.2),
    guide = FALSE
  ) +
  theme(
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    x = "probability of G12D from mutational signatures",
    y = "observed KRAS allele",
    title = "Probabilities of KRAS G12D mutations in COAD tumors samples with various observed KRAS alleles"
  )





not_allele_probs <- all_ranked_allele_predictions %>%
  filter(kras_allele != real_kras_allele) %>%
  group_by(cancer, kras_allele) %>%
  summarise(
    avg_not_allele_prob = mean(allele_prob),
    sd_not_allele_prob = sd(allele_prob)
  ) %>%
  ungroup() %>%
  mutate(
    avg_not_allele_prob_up = avg_not_allele_prob + sd_not_allele_prob,
    avg_not_allele_prob_dn = avg_not_allele_prob - sd_not_allele_prob
  )


is_allele_probs <- all_ranked_allele_predictions %>%
  filter(num_allele_tsb > 10) %>%
  filter(kras_allele == real_kras_allele) %>%
  select(-real_kras_allele) %>%
  group_by(cancer, kras_allele) %>%
  summarise(
    avg_allele_prob = mean(allele_prob),
    sd_allele_prob = sd(allele_prob),
    q25 = quantile(allele_prob, 0.25),
    q75 = quantile(allele_prob, 0.75)
  ) %>%
  ungroup() %>%
  mutate(
    avg_allele_prob_up = avg_allele_prob + sd_allele_prob,
    avg_allele_prob_dn = avg_allele_prob - sd_allele_prob
  )

left_join(is_allele_probs, not_allele_probs,
  by = c("cancer", "kras_allele")
) %>%
  ggplot(aes(x = avg_allele_prob, y = avg_not_allele_prob, color = kras_allele)) +
  facet_wrap(~cancer, nrow = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "grey20", size = 0.6) +
  geom_linerange(
    aes(ymin = avg_not_allele_prob_dn, ymax = avg_not_allele_prob_up),
    alpha = 0.3
  ) +
  geom_linerange(
    aes(xmin = avg_allele_prob_dn, xmax = avg_allele_prob_up),
    alpha = 0.3
  ) +
  geom_point() +
  scale_color_manual(values = short_allele_pal, drop = TRUE) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    x = "average probability of the allele in tumor samples with the allele",
    y = "average probability of the allele in tumor samples without the allele",
    color = "KRAS allele",
    title = "Probability of each KRAS allele in tumor samples with and without the allele"
  )
