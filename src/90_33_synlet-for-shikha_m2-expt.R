# Experimenting on Bayesian models for SL analysis for Shikha.

library(mustashe)
library(ggridges)
library(tidybayes)
library(easystats)
library(rstanarm)
library(ggmcmc)
library(bayesplot)

GRAPHS_DIR <- "90_33_synlet-for-shikha_m2-expt"
reset_graph_directory(GRAPHS_DIR)

theme_set(theme_minimal(base_family = "Arial"))

sample_data_file <- table_path(
  "90_31_synlet-for-shikha_m2",
  "sample-modeling-data.tsv"
)
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
  select(hugo_symbol:dep_map_id, gene_effect:kras_allele,
    is_mutated,
    rna = rna_expression_std
  ) %>%
  mutate(
    kras_allele = factor(kras_allele, levels = kras_allele_order),
    is_mutated = as.numeric(is_mutated)
  )

d %>%
  distinct(dep_map_id, kras_allele) %>%
  count(kras_allele, name = "num_cell_lines")

p <- d %>%
  ggplot(aes(y = hugo_symbol, x = gene_effect)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(size = 0.8, width = 0, height = 0.3, alpha = 0.6)
ggsave_wrapper(p, plot_path(GRAPHS_DIR, "gene-effect_gene.svg"),
  width = 8, height = 10, units = "in"
)

p <- d %>%
  ggplot(aes(x = rna, y = gene_effect)) +
  facet_wrap(~hugo_symbol, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x")
ggsave_wrapper(p, plot_path(GRAPHS_DIR, "rna_gene-effect.svg"),
  width = 12, height = 10, units = "in"
)

p <- d %>%
  ggplot(aes(x = kras_allele, y = gene_effect, color = kras_allele)) +
  facet_wrap(~hugo_symbol, scales = "free") +
  geom_boxplot(size = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.5, size = 1) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave_wrapper(p, plot_path(GRAPHS_DIR, "kras-alleles.svg"),
  width = 14, height = 10, units = "in"
)


m1_prior <- stan_lmer(
  gene_effect ~ 1 +
    (1 | dep_map_id) +
    (1 + rna + kras_allele + is_mutated | hugo_symbol),
  data = d,
  prior = student_t(df = 2, location = 0, scale = 1),
  prior_intercept = normal(location = 0, scale = 0.5),
  prior_aux = cauchy(),
  prior_covariance = decov(),
  prior_PD = FALSE,
  chains = 4,
  adapt_delta = 0.95,
  iter = 4000,
  cores = 4,
  refresh = 2000
)

describe_posterior(m1_prior, effects = "all") %>% as_tibble()

p <- m1_prior %>%
  spread_draws(b[t, g]) %>%
  ungroup() %>%
  filter(str_detect(g, "dep_map_id")) %>%
  ggplot(aes(x = b)) +
  geom_density_ridges(aes(y = g), alpha = 0.3) +
  scale_x_continuous(limits = c(-0.05, 0.05), expand = c(0, 0))
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "m1_cell-line-varying-int_priors.svg"),
  "medium"
)

p <- m1_prior %>%
  spread_draws(b[t, g]) %>%
  ungroup() %>%
  filter(t == "kras_alleleG12D") %>%
  ggplot(aes(x = b)) +
  geom_density_ridges(aes(y = g), alpha = 0.3)
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "m1_krasg12d-varying-slopes_priors.svg"),
  "medium"
)


stash("m1", depends_on = "d", {
  m1 <- stan_lmer(
    gene_effect ~ 1 +
      (1 | dep_map_id) +
      (1 + rna + kras_allele + is_mutated | hugo_symbol),
    data = d,
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    refresh = 1000
  )
})

stash("m1_loo", depends_on = "m1", {
  m1_loo <- loo(m1, k_threshold = 0.7, cores = 4)
})


# Without `is_mutated`.
stash("m2", depends_on = "d", {
  m2 <- stan_lmer(
    gene_effect ~ 1 +
      (1 | dep_map_id) +
      (1 + rna + kras_allele | hugo_symbol),
    data = d,
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.9999,
    iter = 4000,
    refresh = 1000,
    cores = 4
  )
})

stash("m2_loo", depends_on = "m2", {
  m2_loo <- loo(m2, k_threshold = 0.7, cores = 4)
})


# Without varying intercept on cell line.
stash("m3", depends_on = "d", {
  m3 <- stan_lmer(
    gene_effect ~ 1 + (1 + rna + kras_allele | hugo_symbol),
    data = d,
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    refresh = 1000,
    cores = 4
  )
})

stash("m3_loo", depends_on = "m3", {
  m3_loo <- loo(m3, k_threshold = 0.7, cores = 4)
})


# Without varying intercept on cell line and without RNA.
stash("m4", depends_on = "d", {
  m4 <- stan_lmer(
    gene_effect ~ 1 + (1 + kras_allele | hugo_symbol),
    data = d,
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    refresh = 1000
  )
})

stash("m4_loo", depends_on = "m4", {
  m4_loo <- loo(m4, k_threshold = 0.7, cores = 4)
})


# Without varying intercept on cell line and without KRAS allele.
stash("m5", depends_on = "d", {
  m5 <- stan_lmer(
    gene_effect ~ 1 + (1 + rna | hugo_symbol),
    data = d,
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    refresh = 1000
  )
})

stash("m5_loo", depends_on = "m5", {
  m5_loo <- loo(m5, k_threshold = 0.7, cores = 4)
})


# Without without KRAS allele.
stash("m6", depends_on = "d", {
  m6 <- stan_lmer(
    gene_effect ~ 1 + (1 | dep_map_id) + (1 + rna | hugo_symbol),
    data = d,
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    refresh = 1000
  )
})

stash("m6_loo", depends_on = "m6", {
  m6_loo <- loo(m6, k_threshold = 0.7, cores = 2)
})


loo_res <- loo_compare(m1_loo, m2_loo, m3_loo, m4_loo, m5_loo, m6_loo) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(model = fct_inorder(model))
# >   model   elpd_diff   se_diff elpd_loo se_elpd_loo    p_loo se_p_loo     looic se_looic
# > 1    m3   0.0000000 0.0000000 226.3231    19.08046 62.74859 4.933545 -452.6463 38.16092
# > 2    m4  -0.6427467 2.1338088 225.6804    18.68351 57.24245 4.567084 -451.3608 37.36702
# > 3    m2  -1.5359662 0.5265505 224.7872    19.23572 65.57739 5.176505 -449.5743 38.47144
# > 4    m1  -3.2130961 0.9970503 223.1100    19.19069 66.37236 5.179817 -446.2201 38.38139
# > 5    m5 -31.4818113 9.9391388 194.8413    18.55323 25.26455 2.097861 -389.6827 37.10646
# > 6    m6 -32.5976809 9.9878741 193.7255    18.66284 27.42099 2.269134 -387.4509 37.32568


p <- loo_res %>%
  mutate(
    elpd_low = elpd_diff - se_diff,
    elpd_hi = elpd_diff + se_diff
  ) %>%
  ggplot(aes(model, y = elpd_diff)) +
  geom_hline(yintercept = 0) +
  geom_linerange(aes(ymin = elpd_low, ymax = elpd_hi),
    color = "#5c72ed", size = 1.3
  ) +
  geom_line(group = "a", size = 2, alpha = 0.5, color = "#5c72ed") +
  geom_label(aes(label = model),
    size = 4, fill = "#5c72ed", fontface = "bold", family = "Arial",
    color = "white", label.size = 0
  )

ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "model-comparison-loo.svg"),
  "medium"
)

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
p <- m3_draws %>%
  ungroup() %>%
  inner_join(sampled_terms, by = c("term", "group")) %>%
  mutate(model_term = paste0(term, " ", group)) %>%
  ggplot(aes(x = .iteration, y = b, color = factor(.chain))) +
  facet_wrap(~model_term, nrow = 6, scales = "free") +
  geom_line(alpha = 0.6) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "iteration",
    y = "sampled value",
    color = "chain",
    title = "Trace plots from 12 example model terms"
  )
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "m3_example-trace-plots.jpeg"),
  "large"
)
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
p <- m3_draws %>%
  ungroup() %>%
  filter(str_detect(term, "kras_allele")) %>%
  filter(str_detect(group, "hugo_symbol")) %>%
  mutate(
    kras_allele = str_remove(term, "^kras_allele"),
    hugo_symbol = str_remove(group, "^hugo_symbol:")
  ) %>%
  ggplot(aes(x = kras_allele, y = b)) +
  facet_wrap(~hugo_symbol) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(yintercept = c(-0.1, 0.1), lty = 2, color = "grey70") +
  geom_violin(aes(color = kras_allele, fill = kras_allele),
    alpha = 0.5
  ) +
  geom_boxplot(
    color = "black", fill = NA, size = 0.5, width = 0.1,
    outlier.shape = NA
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
ggsave_wrapper(
  p,
  plot_path(GRAPHS_DIR, "m3_post-prob-dist-krasalleles.svg"),
  "medium"
)


m3_post <- d %>%
  modelr::data_grid(hugo_symbol, kras_allele, rna = 0) %>%
  add_predicted_draws(m3)

p <- m3_post %>%
  group_by(hugo_symbol, kras_allele) %>%
  median_hdi(.prediction, .width = c(0.5, 0.75, 0.89, 0.95)) %>%
  ggplot(aes(x = kras_allele)) +
  geom_pointinterval(aes(y = .prediction, ymin = .lower, ymax = .upper),
    alpha = 0.7, color = "grey30"
  ) +
  geom_jitter(aes(y = gene_effect, color = kras_allele),
    data = d, width = 0.3, alpha = 0.7
  ) +
  facet_wrap(~hugo_symbol, scales = "free") +
  scale_color_brewer(palette = "Dark2")
ggsave_wrapper(p, plot_path(GRAPHS_DIR, "m3_post-prediction.svg"), "large")

post_descr <- describe_posterior(m3,
  effects = "random",
  rope_range = c(-0.1, 0.1)
) %>%
  as_tibble()

post_descr %>%
  filter(str_detect(Parameter, "kras_allele")) %>%
  filter(!str_detect(Parameter, "Sigma")) %>%
  mutate(
    info = str_remove_all(Parameter, "^b|\\[|\\]"),
    info = str_remove_all(info, "kras_allele|hugo_symbol|:")
  ) %>%
  separate(info, c("kras_allele", "hugo_symbol"), sep = " ") %>%
  select(kras_allele, hugo_symbol, tidyselect::everything()) %>%
  select(-Parameter) %>%
  arrange(ROPE_Percentage)
