library(glue)
library(magrittr)
library(tidyverse)
library(conflicted)


num_simulations <- 10

num_tumor_samples <- 100

# Get realistic grid of values based on Shamil's paper.
g1_mut_rate <- 0.1
g2_mut_rate <- 0.2


# Get this for actual samples. Maybe a log-transformed rate of mutation or fit
# a gamma distribution.
mean_mutations_per_sample <- 100
sd_mutations_per_sample <- 10


sample_total_mutation_count <- function(n) {
    rnorm(n, mean = mean_mutations_per_sample, sd = sd_mutations_per_sample)
}
# plot(hist(sample_total_mutation_count(1e3)))


sample_gene_mutation <- function(n, mut_rate) {
    rpois(n, mut_rate) > 0
}
# plot(table(sample_gene_mutation(1e3, 0.1)))


tibble(tumor_sample = 1:num_tumor_samples) %>%
    mutate(
        total_num_muts = sample_total_mutation_count(n()),
        g1_mut = sample_gene_mutation(n(), g1_mut_rate),
        g2_mut = sample_gene_mutation(n(), g1_mut_rate))


# calculate OR
# Fisher's test
# count number of False positives
# plot: for various values of OR and number of samples (facet_grid),
#   plot rates of mutation on x and y and heatmap of FDR


# BELOW: experimenting with mutation rates from Shamil's paper

omega <- 0.22
tau <- 01.52
alpha <- 2.23
beta <- 0.37

lambda_s <- seq(0, 2, 0.01)

p <- omega * dexp(lambda_s, tau) + (1-omega) * dinvgamma(lambda_s, alpha, beta)
tibble(p, lambda_s) %>%
    ggplot(aes(x = lambda_s, y = p)) +
    geom_line(group="a")


sample_mutation_rate <- function(n, omega, tau, alpha, beta) {
    omega * rexp(n, tau) + (1-omega) * invgamma::rinvgamma(n, alpha, beta)
}

p_samples <- sample_mutation_rate(1e6, omega, tau, alpha, beta)
tibble(p_samples) %>%
    ggplot(aes(x = p_samples)) +
    geom_histogram(aes(y=..density..), binwidth = 0.01) +
    geom_line(aes(x = lambda_s, y = p), data = tibble(p, lambda_s), group = "a") +
    scale_x_continuous(limits = c(min(lambda_s), max(lambda_s)), expand = c(0, 0))

