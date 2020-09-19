library(mustashe)
library(jhcutils)
library(glue)
library(ggtext)
library(ggridges)
library(patchwork)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
  if (str_detect(f, "global|enrich")) {
    next
  }
  source(f)
}

options(dplyr.summarise.inform = FALSE)

filter_kras_allele_tested <- function(df) {
  df %>%
    filter(is_tested | real_kras_allele == "WT")
}

all_ranked_allele_predictions <- readRDS("~/Downloads/all_ranked_allele_predictions.rds")
ranked_allele_predictions <- readRDS("~/Downloads/ranked_allele_predictions.rds")

theme_set(
  theme_bw(base_size = 11, base_family = "Arial") %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.ticks = element_blank()
    )
)

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"



make_axis_label <- function(final_lbl, ..., sep = "___") {
  addons <- list(...)
  names(addons) <- as.character(seq(1, length(addons)))
  addons <- as_tibble(addons) %>%
    pmap_chr(paste, sep = "-")
  paste(final_lbl, addons, sep = sep)
}


fix_axis_label <- function(x, pattern = "___") {
  unlist(str_split_fixed(x, pattern = pattern, 2)[, 1])
}


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
    x_val = fct_reorder(x_val, accuracy)
  ) %>%
  ggplot(aes(x = x_val, y = accuracy)) +
  facet_grid(~cancer, scales = "free_x", space = "free_x") +
  geom_col() +
  scale_x_discrete(
    expand = c(0, 0),
    labels = fix_axis_label
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  labs(
    x = "observed KRAS allele",
    y = "accuracy"
  )
ggsave_wrapper(
  allele_accuracy_barplots,
  plot_path(GRAPHS_DIR, "allele-accuracy-barplots.svg"),
  "wide"
)










average_allele_probs <- ranked_allele_predictions %>%
  filter_kras_allele_tested() %>%
  group_by(cancer, kras_allele, real_kras_allele) %>%
  summarise(
    avg_allele_prob = mean(allele_prob),
    sd_allele_prob = sd(allele_prob)
  ) %>%
  ungroup()

average_allele_probs_arrows <- average_allele_probs %>%
  select(-sd_allele_prob) %>%
  mutate(
    arrow_alpha = as.numeric(kras_allele == real_kras_allele),
    avg_allele_prob_plus = avg_allele_prob + (max(avg_allele_prob) * 0.1)
  ) %>%
  pivot_longer(c(avg_allele_prob, avg_allele_prob_plus),
               values_to = "y_vals")


pos <- position_dodge(width = 0.9)

average_allele_probs %>%
  mutate(carrot_lbl = "↓",
         carrot_alpha = as.numeric(kras_allele == real_kras_allele)) %>%
  group_by(cancer) %>%
  mutate(carrot_y = avg_allele_prob + (max(avg_allele_prob) * 0.05)) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = kras_allele,
      y = avg_allele_prob,
      fill = real_kras_allele
    )
  ) +
  facet_wrap(~cancer, scales = "free") +
  geom_col(position = pos) +
  geom_text(
    aes(
      label = carrot_lbl,
      y = carrot_y,
      alpha = carrot_alpha
    ),
    position = pos
  ) +
  scale_fill_manual(values = short_allele_pal, drop = TRUE) +
  scale_alpha_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    legend.key.size = unit(3, "mm")
  )
  labs(x = "possible KRAS allele",
       y = "probability",
       fill = "observed KRAS allele",
       color = "observed KRAS allele")
