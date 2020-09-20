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
