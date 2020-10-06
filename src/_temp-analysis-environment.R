library(mustashe)
library(jhcutils)
library(glue)
library(ggridges)
library(patchwork)
library(ggtext)
library(nakedpipe)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
  if (str_detect(f, "global|enrich")) {
    next
  }
  source(f)
}

options(dplyr.summarise.inform = FALSE)

theme_set(theme_bw(base_size = 11, base_family = "Arial"))

ranked_allele_predictions <- readRDS("~/Downloads/ranked_allele_predictions.rds")


################################################################################
# COPIED
GRAPHS_DIR <- "50_13_observed-predicted-kras-alleles_v3_per_tsb"

filter_kras_allele_tested <- function(df) {
  df %>%
    filter(is_tested | real_kras_allele == "WT")
}


################################################################################
