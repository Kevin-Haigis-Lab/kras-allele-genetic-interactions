
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
    theme(axis.ticks = element_blank(),
          legend.position = c(0.8, 0.8),
          legend.background = element_rect(fill = "white", color = "grey50"),
          legend.key.size = unit(3, "mm")) +
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
    plot_path(GRAPHS_DIR,
              glue("allele-prob-other-kras-tumors_{cancer}_{allele}.svg")),
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
            by = c("cancer", "real_kras_allele" = "kras_allele")) %>%
  mutate(adj_allele_prob = allele_prob - avg_incorrect_allele_prob) %>%
  ggplot(aes(x = adj_allele_prob)) +
  facet_wrap(~ cancer, scales = "free", nrow = 2) +
  geom_vline(xintercept = 0, lty = 2, size = 0.6, color = "grey20") +
  geom_density(
    aes(color = real_kras_allele, fill = real_kras_allele),
    alpha = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = short_allele_pal, drop = TRUE) +
  scale_fill_manual(values = short_allele_pal, drop = TRUE) +
  theme(strip.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "probability of observed KRAS allele - average probability in tumors with other KRAS allele",
       y = "density",
       title = "The distribution of the difference in probability of a KRAS allele between\nsamples with the KRAS allele and those with another allele")

ggsave_wrapper(
  diff_between_obs_and_other_plt,
  plot_path(GRAPHS_DIR, "diff-between-obs-and-other-density.svg"),
  "medium"
)
