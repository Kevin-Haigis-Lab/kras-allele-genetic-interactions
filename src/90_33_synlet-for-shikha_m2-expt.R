
library(mustashe)
library(glue)
library(magrittr)

library(tidybayes)
library(easystats)
library(rstanarm)
library(ggmcmc)
library(bayesplot)

library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

GRAPHS_DIR <- "graphs/90_33_synlet-for-shikha_m2-expt"
graphs_path <- function(fname) { file.path(GRAPHS_DIR, fname) }

theme_set(theme_minimal())

sample_data_file <- file.path("~", "Desktop", "sample-modeling-data.tsv")
if (file.exists(sample_data_file)) {
    sample_data <- read_tsv(sample_data_file)
} else {
    stop("Sample data file does not exist.")
}

sample_data %<>%
    filter(cancer == "COAD")

all_kras_alleles <- sort(unique(sample_data$kras_allele))
kras_allele_order <- c("WT", all_kras_alleles[all_kras_alleles != "WT"])


d <- sample_data %>%
    filter(cancer == "COAD") %>%
    select(hugo_symbol:dep_map_id, gene_effect:kras_allele,
           is_mutated, rna = rna_expression_std) %>%
    mutate(kras_allele = factor(kras_allele, levels = kras_allele_order),
           is_mutated = as.numeric(is_mutated))

d %>%
    distinct(dep_map_id, kras_allele) %>%
    count(kras_allele, name = "num_cell_lines")

p <- d %>%
    ggplot(aes(y = hugo_symbol, x = gene_effect)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    geom_jitter(size = 0.8, width = 0, height = 0.3, alpha = 0.6)
ggsave(graphs_path("gene-effect_gene.pdf"), plot = p,
       width = 8, height = 10, units = "in")

p <- d %>%
    ggplot(aes(x = rna, y = gene_effect)) +
    facet_wrap(~ hugo_symbol, scales = "free") +
    geom_point() +
    geom_smooth(method = "lm", formula = "y ~ x")
ggsave(graphs_path("rna_gene-effect.pdf"), plot = p,
       width = 12, height = 10, units = "in")

p <- d %>%
    ggplot(aes(x = kras_allele, y = gene_effect, color = kras_allele)) +
    facet_wrap( ~ hugo_symbol, scales = "free") +
    geom_boxplot(size = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.3, height = 0, alpha = 0.5, size = 1) +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(graphs_path("kras-alleles.pdf"), plot = p,
       width = 14, height = 10, units = "in")


stash("m1", depends_on = "d", {
    m1 <- stan_lmer(
        gene_effect ~ 1 +
            (1|dep_map_id) +
            (1 + rna + kras_allele + is_mutated | hugo_symbol),
        data = d,
        prior = student_t(df = 1, location = 0,
                          scale = 2.5, autoscale = FALSE),
        prior_intercept = normal(location = 0, scale = 1, autoscale = FALSE),
        prior_aux = cauchy(location = 0, scale = 0.5, autoscale = FALSE),
        prior_covariance = decov(regularization = 1, concentration = 1,
                                 shape = 1, scale = 1),
        chains = 4,
        adapt_delta = 0.99,
        iter = 4000,
        cores = 2
    )
})

stash("m1_loo", depends_on = "m1", {
    m1_loo <- loo(m1, k_threshold = 0.7, cores = 2)
})

prior_summary(m1)
#> Priors for model 'm1'
#> ------
#>     Intercept (after predictors centered)
#> ~ normal(location = 0, scale = 1)
#>
#> Auxiliary (sigma)
#> ~ half-cauchy(location = 0, scale = 0.5)
#>
#> Covariance
#> ~ decov(reg. = 1, conc. = 1, shape = 1, scale = 1)
#> ------

# Without `is_mutated`.
stash("m2", depends_on = "d", {
    m2 <- stan_lmer(
        gene_effect ~ 1 +
            (1|dep_map_id) +
            (1 + rna + kras_allele | hugo_symbol),
        data = d,
        prior = student_t(df = 1, location = 0,
                          scale = 2.5, autoscale = FALSE),
        prior_intercept = normal(location = 0, scale = 1, autoscale = FALSE),
        prior_aux = cauchy(location = 0, scale = 0.5, autoscale = FALSE),
        prior_covariance = decov(regularization = 1, concentration = 1,
                                 shape = 1, scale = 1),
        chains = 4,
        adapt_delta = 0.9999,
        iter = 4000,
        cores = 2
    )
})

stash("m2_loo", depends_on = "m2", {
    m2_loo <- loo(m2, k_threshold = 0.7, cores = 2)
})


# Without varying intercept on cell line.
stash("m3", depends_on = "d", {
    m3 <- stan_lmer(
        gene_effect ~ 1 +
            (1 + rna + kras_allele | hugo_symbol),
        data = d,
        prior = student_t(df = 1, location = 0,
                          scale = 2.5, autoscale = FALSE),
        prior_intercept = normal(location = 0, scale = 1, autoscale = FALSE),
        prior_aux = cauchy(location = 0, scale = 0.5, autoscale = FALSE),
        prior_covariance = decov(regularization = 1, concentration = 1,
                                 shape = 1, scale = 1),
        chains = 4,
        adapt_delta = 0.99,
        iter = 4000,
        cores = 2
    )
})

stash("m3_loo", depends_on = "m3", {
    m3_loo <- loo(m3, k_threshold = 0.7, cores = 2)
})


# Without varying intercept on cell line and without RNA.
stash("m4", depends_on = "d", {
    m4 <- stan_lmer(
        gene_effect ~ 1 +
            (1 + kras_allele | hugo_symbol),
        data = d,
        prior = student_t(df = 1, location = 0,
                          scale = 2.5, autoscale = FALSE),
        prior_intercept = normal(location = 0, scale = 1, autoscale = FALSE),
        prior_aux = cauchy(location = 0, scale = 0.5, autoscale = FALSE),
        prior_covariance = decov(regularization = 1, concentration = 1,
                                 shape = 1, scale = 1),
        chains = 4,
        adapt_delta = 0.99,
        iter = 4000,
        cores = 2
    )
})

stash("m4_loo", depends_on = "m4", {
    m4_loo <- loo(m4, k_threshold = 0.7, cores = 2)
})


# Without varying intercept on cell line and without KRAS allele.
stash("m5", depends_on = "d", {
    m5 <- stan_lmer(
        gene_effect ~ 1 +
            (1 + rna | hugo_symbol),
        data = d,
        prior = student_t(df = 1, location = 0,
                          scale = 2.5, autoscale = FALSE),
        prior_intercept = normal(location = 0, scale = 1, autoscale = FALSE),
        prior_aux = cauchy(location = 0, scale = 0.5, autoscale = FALSE),
        prior_covariance = decov(regularization = 1, concentration = 1,
                                 shape = 1, scale = 1),
        chains = 4,
        adapt_delta = 0.99,
        iter = 4000,
        cores = 2
    )
})

stash("m5_loo", depends_on = "m5", {
    m5_loo <- loo(m5, k_threshold = 0.7, cores = 2)
})


# Without without KRAS allele.
stash("m6", depends_on = "d", {
    m6 <- stan_lmer(
        gene_effect ~ 1 +
            (1|dep_map_id) +
            (1 + rna | hugo_symbol),
        data = d,
        prior = student_t(df = 1, location = 0,
                          scale = 2.5, autoscale = FALSE),
        prior_intercept = normal(location = 0, scale = 1, autoscale = FALSE),
        prior_aux = cauchy(location = 0, scale = 0.5, autoscale = FALSE),
        prior_covariance = decov(regularization = 1, concentration = 1,
                                 shape = 1, scale = 1),
        chains = 4,
        adapt_delta = 0.99,
        iter = 4000,
        cores = 2
    )
})

stash("m6_loo", depends_on = "m6", {
    m6_loo <- loo(m6, k_threshold = 0.7, cores = 2)
})


loo_res <- loo_compare(m1_loo, m2_loo, m3_loo, m4_loo, m5_loo, m6_loo) %>%
    as.data.frame() %>%
    rownames_to_column("model") %>%
    mutate(model = fct_inorder(model))
#>   model   elpd_diff   se_diff elpd_loo se_elpd_loo    p_loo se_p_loo     looic se_looic
#> 1    m4   0.0000000  0.000000 226.0296    18.74148 57.33146 4.570858 -452.0591 37.48296
#> 2    m3  -0.3895431  2.116371 225.6400    19.19926 63.60608 5.063226 -451.2800 38.39851
#> 3    m2  -1.4220054  2.214717 224.6075    19.36614 66.03577 5.268828 -449.2151 38.73228
#> 4    m1  -2.3827415  2.198457 223.6468    19.10669 65.86312 5.077927 -447.2936 38.21338
#> 5    m5 -31.4951145  9.963814 194.5344    18.56226 25.54530 2.146073 -389.0689 37.12452
#> 6    m6 -32.2282213 10.009953 193.8013    18.64327 27.30631 2.243467 -387.6027 37.28655

loo_res %>%
    mutate(elpd_low = elpd_diff - se_diff,
           elpd_hi= elpd_diff + se_diff) %>%
    ggplot(aes(model, y = elpd_diff)) +
    geom_hline(yintercept = 0) +
    geom_linerange(aes(ymin = elpd_low, ymax = elpd_hi),
                   color = "#5c72ed", size = 1.3) +
    geom_line(group = "a", size = 2, alpha = 0.5, color = "#5c72ed") +
    geom_label(aes(label = model),
               size = 4, fill = "#5c72ed", fontface = "bold",
               color = "white", label.size = 0)

# launch_shinystan(m1)


head(get_variables(m3))
sample(get_variables(m3), 20)

m3_draws <- m3 %>% spread_draws(b[term, group])

# Example trace plots.
sampled_terms <- m3_draws %>%
    ungroup() %>%
    distinct(term, group) %>%
    group_by(term) %>%
    sample_n(2)
m3_draws %>%
    ungroup() %>%
    inner_join(sampled_terms, by = c("term", "group")) %>%
    mutate(model_term = paste0(term, " ", group)) %>%
    ggplot(aes(x = .iteration, y = b, color = factor(.chain))) +
    facet_wrap(~ model_term, nrow = 6) +
    geom_line(alpha = 0.6) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "iteration",
         y = "sampled value",
         color = "chain",
         title = "Trace plots from 12 example model terms")

# Varying intercept by cell line.
# m3_draws %>%
#     ungroup() %>%
#     filter(term == "(Intercept)" & str_detect(group, "dep_map_id")) %>%
#     mutate(dep_map_id = str_remove(group, "^dep_map_id:")) %>%
#     ggplot(aes(x = b)) +
#     facet_wrap(~ dep_map_id) +
#     geom_vline(xintercept = 0, lty = 2) +
#     geom_density(aes(color = factor(.chain)), size = 1) +
#     scale_color_brewer(palette = "Set1") +
#     labs(x = "posterior values",
#          y = "probability density",
#          color = "chain",
#          title = "Posterior probabilities for varying intercept of cell lines")

# Posterior probabilities for KRAS alleles
m3_draws %>%
    ungroup() %>%
    filter(str_detect(term, "kras_allele")) %>%
    filter(str_detect(group, "hugo_symbol")) %>%
    mutate(kras_allele = str_remove(term, "^kras_allele"),
           hugo_symbol = str_remove(group, "^hugo_symbol:")) %>%
    ggplot(aes(x = kras_allele, y = b)) +
    facet_wrap(~ hugo_symbol) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_hline(yintercept = c(-0.1, 0.1), lty = 2, color = "grey70") +
    geom_violin(aes(color = kras_allele, fill = kras_allele),
                alpha = 0.5) +
    geom_boxplot(color = "black", fill = NA, size = 0.5, width = 0.1,
                 outlier.shape = NA) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")

m3_post <- d %>%
    modelr::data_grid(hugo_symbol, kras_allele, rna = 0) %>%
    add_predicted_draws(m3)

m3_post %>%
    group_by(hugo_symbol, kras_allele) %>%
    median_hdi(.prediction, .width = c(0.5, 0.75, 0.89, 0.95)) %>%
    ggplot(aes(x = kras_allele)) +
    geom_pointinterval(aes(y = .prediction, ymin = .lower, ymax = .upper),
                       alpha = 0.7, color = "grey30") +
    geom_jitter(aes(y = gene_effect, color = kras_allele),
                data = d, width = 0.3, alpha = 0.7) +
    facet_wrap( ~ hugo_symbol, scales = "free") +
    scale_color_brewer(palette = "Dark2")


post_descr <- describe_posterior(m3, effects = "random",
                                 rope_range = c(-0.1, 0.1)) %>%
    as_tibble()

post_descr %>%
    filter(str_detect(Parameter, "kras_allele")) %>%
    filter(!str_detect(Parameter, "Sigma")) %>%
    mutate(info = str_remove_all(Parameter, "^b|\\[|\\]"),
           info = str_remove_all(info, "kras_allele|hugo_symbol|:")) %>%
    separate(info, c("kras_allele", "hugo_symbol"), sep = " ") %>%
    select(kras_allele, hugo_symbol, tidyselect::everything()) %>%
    select(-Parameter) %>%
    arrange(ROPE_Percentage)
