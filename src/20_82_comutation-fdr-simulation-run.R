#!/usr/bin/env Rscript

library(glue)
library(tictoc)
library(magrittr)
library(patchwork)
library(modelr)
library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")

source(file.path("src", "20_83_comutation-fdr-simulation-shared.R"))


sample_gene_mutation <- function(n, mut_rate) {
    rpois(n, mut_rate) > 0
}

# Ensure there is always a 2x2 contingency table.
my_table <- function(a, b) {
    a <- factor(a, levels = c("FALSE", "TRUE"))
    b <- factor(b, levels = c("FALSE", "TRUE"))
    table(a, b)
}

comutation_simulation <- function(num_tumor_samples, g1_rate, g2_rate, ...) {
    mut_data <- tibble(tumor_sample = seq(1, num_tumor_samples)) %>%
        mutate(
            g1_mut = sample_gene_mutation(n(), g1_rate),
            g2_mut = sample_gene_mutation(n(), g2_rate))

    mut_table <- my_table(mut_data$g1_mut, mut_data$g2_mut)
    fish_res <- fisher.test(mut_table, alternative = "g") %>%
        broom::tidy()
    return(bind_cols(
        tibble(mutation_table = list(mut_table)),
        fish_res
    ))
}

args <- commandArgs(trailingOnly = TRUE)
job_array_idx <- args[[1]]

message(glue("Running simulation job index {job_array_idx}."))

tic("simulation")
simulation_df <- read_tsv(simulation_input_file_name(job_array_idx))
simulation_df$simulation_res <- pmap(simulation_df, comutation_simulation)
simulation_df <- simulation_df %>%
    unnest(simulation_res) %>%
    janitor::clean_names()
toc()

qs::qsave(simulation_df, simulation_output_file_name(job_array_idx))

message("Finished simulation.")
