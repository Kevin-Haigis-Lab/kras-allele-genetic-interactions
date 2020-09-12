
library(mustashe)
library(jhcutils)
library(glue)
library(ggtext)
library(patchwork)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
  if (str_detect(f, "global|enrich")) {
    next
  }
  source(f)
}

options(dplyr.summarise.inform = FALSE)

real_kras_mutations <- readRDS("~/Downloads/real_kras_mutations.rds")
kras_allele_predictions <- readRDS("~/Downloads/kras_allele_predictions.rds")
all_kras_allele_predictions <- readRDS("~/Downloads/all_kras_allele_predictions.rds")


# Filter out infrequent KRAS alleles
filter_kras_allele_ct <- function(df, min_n = 20) {
    original_groups <- dplyr::groups(df)
    df %>%
        group_by(cancer, real_kras_allele) %>%
        filter(n_distinct(tumor_sample_barcode) >= !!min_n) %>%
        group_by(!!!original_groups)
}


ranked_allele_predictions <- kras_allele_predictions %>%
  inner_join(
    real_kras_mutations %>% rename(real_kras_allele = kras_allele),
    by = c("cancer", "tumor_sample_barcode")
  ) %>%
  group_by(cancer, tumor_sample_barcode) %>%
  arrange(-allele_prob) %>%
  mutate(allele_idx = row_number()) %>%
  ungroup() %>%
  arrange(cancer, tumor_sample_barcode, allele_idx)

theme_set(theme_bw(base_size = 7, base_family = "Arial"))

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"




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
  title <- glue("Probability of a {allele} allele in {cancer} tumors with various KRAS alleles")

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
      legend.position = c(0.8, 0.8),
      legend.background = element_rect(fill = "white", color = "grey50"),
      legend.key.size = unit(3, "mm")
    ) +
    labs(
      x = xlbl,
      y = "density",
      color = "observed KRAS allele\nof tumor",
      fill = "observed KRAS allele\nof tumor",
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



calculate_frac_with_lower_prob <- function(prob, non_mut_df) {
  mean(non_mut_df$allele_prob < prob)
}


allele_prob_dist_plot1 <- function(df, allele) {
  pal <- c("TRUE" = short_allele_pal[[allele]],
           "FALSE" = "grey50")
  lbls <- c("TRUE" = glue("{allele} mutant"),
            "FALSE" = glue("non-{allele} mutant"))
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
    theme(axis.ticks = element_blank(),
          legend.position = c(0.75, 0.8),
          legend.key.size = unit(3, "mm")) +
    labs(x = xlbl,
         y = "density",
         color = NULL,
         fill = NULL)
}


allele_prob_frac_dist_plot2 <- function(df, allele) {
  xlbl <- glue("fraction of non-{allele} mutant samples\nwith a lower probability of {allele} mutation")
  df %>%
    ggplot(aes(x = frac_less_prob)) +
    geom_vline(xintercept = 0.5, lty = 2, color = "grey20", size = 0.6) +
    geom_density(fill = "#8491a1", alpha = 0.5, size = 1, color = "#8491a1") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = xlbl,
         y = "density")
}


fraction_of_allele_prob_plot <- function(cancer,
                                         allele,
                                         prob_data,
                                         ignore_alleles = NULL) {
  prob_allele <- prob_data %>%
    filter(cancer == !!cancer) %>%
    filter(kras_allele == !!allele & real_kras_allele == !!allele)

  prob_not_allele <- prob_data %>%
    filter(cancer == !!cancer) %>%
    filter(kras_allele == !!allele & real_kras_allele != !!allele) %>%
    filter(!(kras_allele %in% !!ignore_alleles))

  p1 <- prob_data %>%
    filter(cancer == !!cancer & kras_allele == !!allele) %>%
    mutate(is_allele = real_kras_allele == !!allele) %>%
    allele_prob_dist_plot1(allele = allele)

  p2 <- prob_allele %>%
    mutate(frac_less_prob = map_dbl(allele_prob,
                                    calculate_frac_with_lower_prob,
                                    non_mut_df = prob_not_allele)) %>%
    allele_prob_frac_dist_plot2(allele = allele)

  patch_title <- glue("The probability of {allele} mutations in {cancer} tumors with and without a {allele} allele")
  patch <- (p1 | p2) +
    plot_annotation(title = patch_title,
                    theme = theme(plot.title = element_text(hjust = 0.5)))

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
  left_join(tibble(cancer = "LUAD",
                   allele = c("G12C", "G13C"),
                   ignore_alleles = c("G13C", "G12C")),
            by = c("cancer", "allele")) %>%
  pwalk(fraction_of_allele_prob_plot, prob_data = ranked_allele_predictions)
