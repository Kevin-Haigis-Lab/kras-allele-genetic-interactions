# Analysis of the predicted KRAS alleles per sample.

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"
reset_graph_directory(GRAPHS_DIR)

library(ggridges)

set.seed(0)

theme_set(theme_bw(base_size = 7, base_family = "Arial"))

# Minumum number of tumor samples to include the KRAS allele in the analysis.
MIN_ALLELE_TSB <- 15

filter_kras_allele_ct <- function(df) {
  df %>%
    filter(num_allele_tsb >= MIN_ALLELE_TSB)
}

filter_kras_allele_tested <- function(df) {
  df %>%
    filter(is_tested | real_kras_allele == "WT")
}

add_kras_allele_tested <- function(df) {
  df %>%
    left_join(
      alleles_for_each_cancer_obs_vs_pred %>%
        add_column(is_tested = TRUE),
      by = c("cancer", "real_kras_allele" = "kras_allele")
    ) %>%
    mutate(is_tested = ifelse(is.na(is_tested), FALSE, is_tested))
}


real_kras_mutations <- trinucleotide_mutations_df %>%
  distinct(cancer, tumor_sample_barcode, kras_allele)


prepare_ranked_allele_prediction_df <- function(df) {
  df %>%
    inner_join(
      real_kras_mutations %>%
        rename(real_kras_allele = kras_allele),
      by = c("tumor_sample_barcode", "cancer")
    ) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    arrange(-allele_prob) %>%
    mutate(allele_idx = row_number()) %>%
    group_by(cancer, real_kras_allele) %>%
    mutate(num_allele_tsb = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    arrange(cancer, tumor_sample_barcode, allele_idx) %>%
    add_kras_allele_tested()
}


cache(
  "ranked_allele_predictions",
  depends = "kras_allele_predictions",
  {
    ranked_allele_predictions <- prepare_ranked_allele_prediction_df(
      kras_allele_predictions
    )
    return(ranked_allele_predictions)
  }
)


cache(
  "all_ranked_allele_predictions",
  depends = "all_kras_allele_predictions",
  {
    all_ranked_allele_predictions <- prepare_ranked_allele_prediction_df(
      all_kras_allele_predictions
    )
    return(all_ranked_allele_predictions)
  }
)



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
    fill_idx_lbls <- rev(as.character(seq(1, max(df$allele_idx))))
  }

  df %>%
    filter(real_kras_allele == kras_allele) %>%
    mutate(allele_idx = minmax(allele_idx, 0, max_pred_idx)) %>%
    count(cancer, real_kras_allele, num_allele_tsb, allele_idx) %>%
    mutate(
      correct_allele = allele_idx == 1,
      real_kras_allele = factor_alleles(real_kras_allele)
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
      title = "The rank of the observed allele when ordered by probability using mutational signatures"
    )
}

ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  plot_allele_probability_barplot() %>%
  ggsave_wrapper(
    plot_path(GRAPHS_DIR, "allele-prob-ranking_stack_barplot.svg"),
    "medium"
  )

p <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  plot_allele_probability_barplot(geombar_position = "fill")
p <- p +
  scale_y_continuous(expand = c(0, 0))
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "allele-prob-ranking_fill_barplot.svg"),
  "medium"
)


all_ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  plot_allele_probability_barplot(max_pred_idx = 7) %>%
  ggsave_wrapper(
    plot_path(GRAPHS_DIR, "all-allele-prob-ranking_stack_barplot.svg"),
    "medium"
  )

p <- all_ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  plot_allele_probability_barplot(geombar_position = "fill", max_pred_idx = 7)
p <- p +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "fraction of tumor samples")
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
  fname <- glue("ranked-probability-lines_{cancer}_{allele}.svg")
  ggsave_wrapper(
    plt,
    plot_path(GRAPHS_DIR, fname),
    "wide"
  )
}

ranked_prob_plot_data <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
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
  filter_kras_allele_tested() %>%
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

allele_prediction_barplot <- function(df,
                                      y_values,
                                      y_lbl = NULL,
                                      fill_lbl = NULL,
                                      bar_position = "fill") {
  df %>%
    ggplot(aes(x = real_kras_allele, y = {{ y_values }})) +
    facet_wrap(
      ~cancer,
      scales = "free",
      nrow = 2
    ) +
    geom_col(
      aes(fill = kras_allele),
      position = "fill",
      width = 0.8
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = short_allele_pal) +
    theme(
      legend.key.size = unit(3, "mm"),
      plot.title = element_text(hjust = 0.5),
      strip.background = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      x = "observed KRAS allele",
      y = y_lbl,
      fill = fill_lbl
    )
}

walk2(
  list(ranked_allele_predictions, all_ranked_allele_predictions),
  c("most-likely-other-allele.svg", "most-likely-other-allele_all-alleles.svg"),
  function(df, fn) {
    p <- df %>%
      filter_kras_allele_tested() %>%
      filter(allele_idx == 1) %>%
      count(cancer, real_kras_allele, kras_allele) %>%
      mutate(
        real_kras_allele = factor_alleles(real_kras_allele),
        kras_allele = factor_alleles(kras_allele)
      ) %>%
      allele_prediction_barplot(
        y_values = n,
        y_lbl = "fraction of tumor samples most most likely to be another allele",
        fill_lbl = "most likely\nKRAS allele"
      )
    ggsave_wrapper(
      p,
      plot_path(GRAPHS_DIR, fn),
      "wide"
    )
  }
)

walk2(
  list(ranked_allele_predictions, all_ranked_allele_predictions),
  c("cum-prob-of-other-allele.svg", "cum-prob-of-other-allele_all-alleles.svg"),
  function(df, fn) {
    p <- df %>%
      filter_kras_allele_tested() %>%
      group_by(cancer, real_kras_allele, kras_allele) %>%
      summarise(total_prob = sum(allele_prob)) %>%
      ungroup() %>%
      mutate(
        real_kras_allele = factor_alleles(real_kras_allele),
        kras_allele = factor_alleles(kras_allele)
      ) %>%
      allele_prediction_barplot(
        y_values = total_prob,
        y_lbl = "fraction of tumor samples most most likely to be another allele",
        fill_lbl = "alternative\nKRAS allele"
      )
    ggsave_wrapper(
      p,
      plot_path(GRAPHS_DIR, fn),
      "wide"
    )
  }
)



#### ---- Distributions of allele probabilities ---- ####

# All in one large, faceted plot.
distribution_of_allele_probs <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
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


# Individual ggridge plots per allele per cancer.

allele_prob_distribution_per_kras_mutation <- function(df, cancer, allele) {
  xlbl <- glue("probability of {allele} allele")
  title <- glue("Probability of a {allele} allele in {cancer} tumors\nwith various KRAS alleles")

  df %>%
    filter_kras_allele_tested() %>%
    filter(cancer == !!cancer & kras_allele == !!allele) %>%
    mutate(is_allele = real_kras_allele == !!allele) %>%
    ggplot(aes(x = allele_prob, y = real_kras_allele)) +
    geom_density_ridges(
      aes(
        color = real_kras_allele,
        fill = real_kras_allele,
        alpha = is_allele
      ),
      adjust = 1.5
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
      values = c("TRUE" = 0.9, "FALSE" = 0.3),
      guide = FALSE
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_discrete(expand = expansion(mult = c(0, 0.02))) +
    theme(
      axis.ticks = element_blank()
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
  filter_kras_allele_tested() %>%
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


#### ---- Ratio of allele prob to average ---- ####

average_allele_prob_incorrect_alleles <- ranked_allele_predictions %>%
  filter(kras_allele != real_kras_allele) %>%
  group_by(cancer, kras_allele) %>%
  summarise(avg_incorrect_allele_prob = mean(allele_prob)) %>%
  ungroup()

diff_between_obs_and_other_plt <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
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


#### ---- Fraction with lower probability in tumors with other alleles ---- ####


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



#### ---- Plots for figures ---- ####


codon_pal <- codon_palette[names(codon_palette) != "Other"]

ranked_allele_predictions_top <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  filter(allele_idx == 1) %>%
  mutate(is_correct = kras_allele == real_kras_allele)

allele_predictions_acc <- ranked_allele_predictions_top %>%
  filter(real_kras_allele != "WT") %>%
  group_by(cancer, real_kras_allele) %>%
  summarise(accuracy = mean(as.numeric(is_correct))) %>%
  ungroup() %>%
  mutate(
    false_pos_mut = map2_dbl(
      cancer,
      real_kras_allele,
      function(C, A) {
        ranked_allele_predictions_top %>%
          filter(cancer == C & real_kras_allele != A) %>%
          filter(real_kras_allele != "WT") %>%
          mutate(is_wrong = kras_allele == A) %>%
          pull(is_wrong) %>%
          unlist() %>%
          mean()
      }
    ),
    false_pos_wt = map2_dbl(
      cancer,
      real_kras_allele,
      function(C, A) {
        ranked_allele_predictions_top %>%
          filter(cancer == C & real_kras_allele == "WT") %>%
          mutate(is_wrong = kras_allele == A) %>%
          pull(is_wrong) %>%
          unlist() %>%
          mean()
      }
    )
  )

knitr::kable(allele_predictions_acc, format = "markdown", digits = 2)
# > |cancer |real_kras_allele | accuracy| false_pos_mut| false_pos_wt|
# > |:------|:----------------|--------:|-------------:|------------:|
# > |COAD   |A146T            |     0.17|          0.20|         0.16|
# > |COAD   |G12A             |     0.05|          0.04|         0.03|
# > |COAD   |G12C             |     0.08|          0.05|         0.03|
# > |COAD   |G12D             |     0.17|          0.22|         0.19|
# > |COAD   |G12S             |     0.00|          0.08|         0.08|
# > |COAD   |G12V             |     0.14|          0.07|         0.07|
# > |COAD   |G13D             |     0.33|          0.34|         0.43|
# > |LUAD   |G12A             |     0.00|          0.01|         0.04|
# > |LUAD   |G12C             |     0.62|          0.61|         0.47|
# > |LUAD   |G12D             |     0.04|          0.05|         0.19|
# > |LUAD   |G12V             |     0.26|          0.35|         0.30|
# > |LUAD   |G13C             |     0.00|          0.00|         0.00|
# > |MM     |G12A             |     0.00|          0.05|         0.08|
# > |MM     |G12D             |     0.38|          0.24|         0.26|
# > |MM     |G12R             |     0.00|          0.00|         0.01|
# > |MM     |G12V             |     0.16|          0.08|         0.09|
# > |MM     |G13D             |     0.38|          0.34|         0.33|
# > |MM     |Q61H             |     0.20|          0.21|         0.18|
# > |MM     |Q61L             |     0.17|          0.02|         0.01|
# > |MM     |Q61R             |     0.00|          0.02|         0.04|
# > |PAAD   |G12C             |     0.15|          0.08|         0.07|
# > |PAAD   |G12D             |     0.32|          0.28|         0.30|
# > |PAAD   |G12R             |     0.00|          0.01|         0.03|
# > |PAAD   |G12V             |     0.39|          0.40|         0.32|
# > |PAAD   |Q61H             |     0.20|          0.22|         0.28|

allele_accuracy_barplots <- allele_predictions_acc %>%
  mutate(
    x_val = factor_alleles(real_kras_allele),
    codon = str_extract(real_kras_allele, "[:digit:]+"),
    codon = fct_reorder(codon, as.numeric(codon), .fun = unique)
  ) %>%
  ggplot(aes(x = x_val, y = accuracy)) +
  facet_grid(. ~ cancer, space = "free_x", scales = "free_x") +
  geom_linerange(
    aes(ymin = 0, ymax = accuracy),
    size = 1,
    color = "grey75"
  ) +
  geom_point(
    aes(color = codon),
    size = 2
  ) +
  geom_point(
    aes(y = false_pos_mut, alpha = cancer),
    size = 1.3,
    shape = 21,
    color = "grey35"
  ) +
  geom_point(
    aes(y = false_pos_wt),
    size = 1.3,
    shape = 25,
    color = "grey35"
  ) +
  scale_y_continuous(
    limits = c(0, 0.7),
    breaks = seq(0, 0.7, 0.1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_color_manual(
    values = codon_pal,
    guide = guide_legend(
      order = 10,
      override.aes = list(size = 1.3, alpha = 1)
    )
  ) +
  scale_alpha_manual(
    values = c(1, 1, 1, 1),
    breaks = c("COAD", "LUAD", "MM"),
    labels = c("the *KRAS* allele", "another *KRAS* mutation", "*KRAS* WT"),
    guide = guide_legend(
      title = "tumor samples with",
      override.aes = list(
        shape = c(21, 21, 25),
        size = c(2, 1.3, 1.3),
        alpha = 1,
        color = "grey20",
        fill = c("grey20", "white", "white")
      )
    )
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_markdown()
  ) +
  labs(
    x = NULL,
    y = "fraction of tumor samples\npredicted to have allele"
  )
ggsave_wrapper(
  allele_accuracy_barplots,
  plot_path(GRAPHS_DIR, "allele-accuracy-barplots.svg"),
  "wide"
)
saveFigRds(allele_accuracy_barplots, "allele_accuracy_barplots")




average_allele_probs <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  group_by(cancer, kras_allele, real_kras_allele) %>%
  summarise(
    avg_allele_prob = median(allele_prob),
    sd_allele_prob = sd(allele_prob),
    allele_prob_q25 = quantile(allele_prob, 0.25),
    allele_prob_q75 = quantile(allele_prob, 0.75)
  ) %>%
  ungroup()

average_allele_lines <- ranked_allele_predictions %>%
  mutate(
    is_allele = real_kras_allele == kras_allele,
    is_allele = fct_rev(factor(is_allele)),
    kras_allele = factor_alleles(kras_allele)
  ) %>%
  group_by(cancer, kras_allele, is_allele) %>%
  summarise(avg_allele_prob = median(allele_prob)) %>%
  ungroup()

wide_average_allele_lines <- average_allele_lines %>%
  mutate(is_allele = ifelse(is_allele == "TRUE", "ymax", "ymin")) %>%
  pivot_wider(
    c(cancer, kras_allele),
    names_from = is_allele,
    values_from = avg_allele_prob
  )

pos <- position_dodge(width = 0.7)

allele_prob_barplot_arrows <- average_allele_probs %>%
  mutate(
    carrot_lbl = "↓",
    carrot_alpha = as.numeric(kras_allele == real_kras_allele),
    y_up = avg_allele_prob + sd_allele_prob,
    y_dn = avg_allele_prob - sd_allele_prob,
    y_dn = minmax(y_dn, 0, 100),
    real_kras_allele = factor_alleles(real_kras_allele),
    kras_allele = factor_alleles(kras_allele)
  ) %>%
  group_by(cancer) %>%
  mutate(carrot_y = avg_allele_prob + (max(avg_allele_prob) * 0.1)) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = kras_allele,
      y = avg_allele_prob,
      color = real_kras_allele
    )
  ) +
  facet_wrap(~cancer, nrow = 2, scales = "free") +
  geom_crossbar(
    aes(ymin = ymin, ymax = ymax, y = ymin),
    fill = "grey50",
    color = NA,
    alpha = 0.2,
    data = wide_average_allele_lines
  ) +
  geom_errorbar(
    aes(ymin = avg_allele_prob, ymax = avg_allele_prob, linetype = is_allele),
    data = average_allele_lines,
    color = "grey30",
    alpha = 0.5,
    size = 0.3
  ) +
  geom_linerange(
    aes(
      ymin = allele_prob_q25,
      ymax = allele_prob_q75
    ),
    alpha = 0.4,
    size = 0.3,
    position = pos
  ) +
  geom_point(
    position = pos,
    alpha = 1,
    size = 0.8
  ) +
  scale_color_manual(
    values = short_allele_pal,
    drop = TRUE,
    guide = guide_legend(
      order = 10,
      override.aes = list(size = 1.3, lty = 0, alpha = 1)
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_linetype_manual(
    values = c("TRUE" = 1, "FALSE" = 6),
    labels = c(
      "TRUE" = "with the\nKRAS allele",
      "FALSE" = "with another\nKRAS allele"
    ),
    guide = guide_legend(
      override.aes = list(
        lty = c("TRUE" = 1, "FALSE" = 3),
        alpha = 1,
        size = 0.6
      ),
      keyheight = unit(6, "mm"),
      order = 20
    )
  ) +
  theme(
    legend.key.size = unit(3, "mm"),
    strip.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "possible KRAS allele",
    y = "probability of KRAS allele",
    color = "observed\nKRAS allele",
    linetype = "average of\ntumor samples"
  )
ggsave_wrapper(
  allele_prob_barplot_arrows,
  plot_path(GRAPHS_DIR, "allele-prob-barplot_arrows.svg"),
  "wide"
)
saveFigRds(allele_prob_barplot_arrows, "allele_prob_barplot_arrows")




is_real_above_other <- function(is_allele, val) {
  stopifnot(length(is_allele) == 2 & length(val) == 2)
  y <- val[order(is_allele)]
  return(y[[1]] < y[[2]])
}


ifelse_pal <- function(x, options) {
  ifelse(x, names(options)[[1]], names(options)[[2]])
}


boot_wrapper <- function(x, fxn, indices) {
  fxn(x[indices])
}

boot_median <- function(x, indices) {
  boot_wrapper(x, median, indices)
}

boot_mean <- function(x, indices) {
  boot_wrapper(x, mean, indices)
}

point_pal <- c(
  "the *KRAS* allele" = "darkslateblue",
  "another *KRAS* mutation" = "grey20",
  "*KRAS* WT" = "grey50"
)

shape_pal <- c(19, 21, 25)
names(shape_pal) <- names(point_pal)


allele_prob_per_allele_df <- ranked_allele_predictions %.% {
  filter_kras_allele_tested()
  mutate(
    allele_group = case_when(
      kras_allele == real_kras_allele ~ "the *KRAS* allele",
      real_kras_allele == "WT" ~ "*KRAS* WT",
      TRUE ~ "another *KRAS* mutation"
    ),
    allele_group = factor(allele_group, levels = names(point_pal)),
    kras_allele = factor_alleles(kras_allele)
  )
}

allele_prob_stats <- function(df) {
  pairwise.wilcox.test(
    x = df$allele_prob,
    g = df$allele_group,
    p.adjust.method = "BH"
  ) %>%
    broom::tidy() %>%
    rename(adj_p_value = p.value)
}

allele_prob_per_allele_stats <- allele_prob_per_allele_df %.% {
  group_by(cancer, kras_allele)
  nest()
  ungroup()
  mutate(stats_res = map(data, allele_prob_stats))
  select(-data)
  unnest(stats_res)
  filter(group1 == "the *KRAS* allele" | group2 == "the *KRAS* allele")
  mutate(
    comparison = ifelse(group1 == "the *KRAS* allele", group2, group1),
    comparison = factor(comparison, levels = names(point_pal))
  )
}

allele_prob_per_allele_summary <- allele_prob_per_allele_df %.% {
  group_by(cancer, allele_group, kras_allele)
  summarise(
    mean_prob = mean(allele_prob),
    median_prob = median(allele_prob),
    mean_prob_boot = list(boot::boot(allele_prob, boot_mean, R = 1e4)),
    q25_prob = quantile(allele_prob, 0.25),
    q75_prob = quantile(allele_prob, 0.75)
  )
  ungroup()
  mutate(
    mean_ci = map(
      mean_prob_boot,
      boot::boot.ci,
      conf = 0.95,
      type = "perc"
    ),
    mean_ci_lower = map_dbl(mean_ci, ~ .x$perc[[4]]),
    mean_ci_upper = map_dbl(mean_ci, ~ .x$perc[[5]])
  )
  left_join(
    allele_prob_per_allele_stats %>% filter(adj_p_value < 0.05),
    by = c("cancer", "kras_allele", "allele_group" = "comparison")
  )
  mutate(star_point_size = !is.na(adj_p_value))
}


pos_width <- 0.7
pos <- position_dodge(width = pos_width)

allele_prob_per_allele_plot_annotations <- allele_prob_per_allele_stats %.% {
  filter(adj_p_value < 0.05)
  left_join(
    allele_prob_per_allele_summary %>%
      filter(allele_group == "the *KRAS* allele") %>%
      select(cancer, kras_allele, mean_ci_upper),
    by = c("cancer", "kras_allele")
  )
  arrange(cancer, kras_allele, group1)
  mutate(
    p_value_label = format_pvalue_label(adj_p_value, prefix = ""),
    p_x = c(5, 5, 2),
    p_y = mean_ci_upper + c(0.05, 0.02, 0.02),
    seg_x1 = p_x - 0.3 * pos_width,
    seg_x2 = p_x + c(0.43, 0, 0.43) * pos_width,
    seg_y2 = p_y - 0.007,
    seg_y1 = seg_y2 - 0.004
  )
}

segment_size <- 0.3

allele_prob_per_allele_plot <- allele_prob_per_allele_summary %>%
  ggplot(aes(x = kras_allele, y = mean_prob)) +
  facet_grid(. ~ cancer, scales = "free", space = "free_x") +
  geom_linerange(
    aes(ymin = mean_ci_lower, ymax = mean_ci_upper, color = allele_group),
    position = pos,
    size = 0.3,
    alpha = 0.75
  ) +
  geom_point(
    aes(color = allele_group, shape = allele_group),
    size = 1.3,
    position = pos,
    fill = "white"
  ) +
  geom_text(
    data = allele_prob_per_allele_plot_annotations,
    aes(x = p_x, y = p_y, label = p_value_label),
    size = 1.6,
    family = "Arial",
    hjust = 0.5,
    vjust = 0
  ) +
  geom_segment(
    data = allele_prob_per_allele_plot_annotations,
    aes(x = seg_x2, xend = seg_x2, y = seg_y1, yend = seg_y2),
    size = segment_size,
    color = "black",
    lineend = "round"
  ) +
  geom_segment(
    data = allele_prob_per_allele_plot_annotations,
    aes(x = seg_x1, xend = seg_x1, y = seg_y1, yend = seg_y2),
    size = segment_size,
    color = "black",
    lineend = "round"
  ) +
  geom_segment(
    data = allele_prob_per_allele_plot_annotations,
    aes(x = seg_x1, xend = seg_x2, y = seg_y2, yend = seg_y2),
    size = segment_size,
    color = "black",
    lineend = "round"
  ) +
  scale_color_manual(
    values = point_pal,
    guide = guide_legend(
      override.aes = list(shape = shape_pal),
      order = 10
    )
  ) +
  scale_shape_manual(
    values = shape_pal,
    guide = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_markdown()
  ) +
  labs(
    x = NULL,
    y = "average probability",
    color = "tumor samples with"
  )
ggsave_wrapper(
  allele_prob_per_allele_plot,
  plot_path(GRAPHS_DIR, "allele-prob-per-allele_scatter.svg"),
  "wide"
)
saveFigRds(allele_prob_per_allele_plot, "allele_prob_per_allele_plot")


stats_stars_plot <- allele_prob_per_allele_stats %.%
  {
    mutate(
      star_lbl = assign_stars(adj_p_value),
      comparison = ifelse(group1 == "the KRAS allele", group2, group1),
      comparison = ifelse(comparison == "another KRAS mutation", "other KRAS mutants", comparison),
      comparison = glue("vs. {comparison}")
    )
  } %>%
  ggplot(aes(kras_allele, comparison)) +
  facet_grid(. ~ cancer, scales = "free", space = "free_x") +
  geom_tile(fill = "white", color = NA) +
  geom_text(
    aes(label = star_lbl),
    family = "Arial",
    size = 2,
    hjust = 0.5,
    vjust = 0.75
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_markdown(),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    axis.line = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  )


allele_prob_per_allele_plot_MOD <- allele_prob_per_allele_plot +
  theme(strip.text = element_blank())

allele_prob_per_allele_patch <- (
  stats_stars_plot / allele_prob_per_allele_plot_MOD
) +
  plot_layout(heights = c(1, 15))

ggsave_wrapper(
  allele_prob_per_allele_patch,
  plot_path(GRAPHS_DIR, "allele-prob-per-allele_with-stats.svg"),
  "wide"
)
