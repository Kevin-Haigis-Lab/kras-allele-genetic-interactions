
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

theme_set(theme_bw(base_size = 11, base_family = "Arial"))

GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"



