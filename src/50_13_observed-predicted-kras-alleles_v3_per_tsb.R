# Analysis of the predicted KRAS alleles per sample.

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"
reset_graph_directory(GRAPHS_DIR)

library(ggridges)

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



allele_predictions_acc <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  filter(real_kras_allele != "WT") %>%
  filter(allele_idx == 1) %>%
  mutate(is_correct = kras_allele == real_kras_allele) %>%
  group_by(cancer, real_kras_allele) %>%
  summarise(accuracy = mean(as.numeric(is_correct))) %>%
  ungroup()

allele_accuracy_barplots <- allele_predictions_acc %>%
  mutate(
    x_val = make_axis_label(real_kras_allele, cancer),
    x_val = fct_reorder(x_val, accuracy),
    codon = str_extract(real_kras_allele, "[:digit:]+"),
    codon = fct_reorder(codon, as.numeric(codon), .fun = unique)
  ) %>%
  ggplot(aes(x = x_val, y = accuracy)) +
  facet_wrap(~cancer, nrow = 1, scales = "free_x") +
  geom_linerange(
    aes(ymin = 0, ymax = accuracy),
    size = 1,
    color = "grey75"
  ) +
  geom_point(
    aes(color = codon),
    size = 2
  ) +
  scale_x_discrete(
    labels = fix_axis_label
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
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "grey30")
  ) +
  labs(
    x = "observed KRAS allele",
    y = "fraction of tumor samples\nwith correctly predicted allele"
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
