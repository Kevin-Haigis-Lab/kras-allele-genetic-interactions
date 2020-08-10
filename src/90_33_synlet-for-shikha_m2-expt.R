
library(mustashe)
library(glue)

library(tidybayes)
library(easystats)
library(rstanarm)
library(ggmcmc)
library(bayesplot)

library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

theme_set(theme_minimal())

sample_data_file <- file.path("~", "Desktop", "sample-modeling-data.tsv")
if (file.exists(sample_data_file)) {
    sample_data <- read_tsv(sample_data_file)
} else {
    stop("Sample data file does not exist.")
}

all_kras_alleles <- sort(unique(sample_data$kras_allele))
kras_allele_order <- c("WT", all_kras_alleles[all_kras_alleles != "WT"])


d <- sample_data %>%
    select(hugo_symbol:dep_map_id, gene_effect:kras_allele,
           is_mutated, rna = rna_expression_std) %>%
    mutate(kras_allele = factor(kras_allele, levels = kras_allele_order),
           is_mutated = as.numeric(is_mutated))

d %>%
    ggplot(aes(y = hugo_symbol, x = gene_effect)) +
    geom_boxplot(aes(color = cancer, fill = cancer),
                 alpha = 0.2, outlier.shape = NA) +
    geom_jitter(aes(color = cancer),
                size = 0.8, width = 0, height = 0.3, alpha = 0.6) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")
ggsave("~/Desktop/gene-effect_gene.pdf", width = 8, height = 10, units = "in")

d %>%
    ggplot(aes(x = rna, y = gene_effect)) +
    facet_wrap(~ hugo_symbol, scales = "free") +
    geom_point(aes(color = cancer)) +
    geom_smooth(method = "lm", formula = "y ~ x") +
    scale_color_brewer(palette = "Set1")
ggsave("~/Desktop/rna_gene-effect.pdf", width = 12, height = 10, units = "in")

d %>%
    ggplot(aes(x = kras_allele, y = gene_effect, color = kras_allele)) +
    facet_wrap(cancer ~ hugo_symbol, scales = "free", nrow = 4) +
    geom_boxplot(size = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.3, height = 0, alpha = 0.5, size = 1) +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("~/Desktop/kras-alleles.pdf", width = 20, height = 10, units = "in")


stash("m1", depends_on = "d", {
    m1 <- stan_lmer(
        gene_effect ~ 1 +
            (1|dep_map_id) +
            (1 + rna + kras_allele + is_mutated | cancer/hugo_symbol),
        data = d,
        cores = 2
    )
})

# launch_shinystan(m1)


get_variables(m1)

m1_draws <- m1 %>% spread_draws(b[term, group])
m1_draws %>%
    ungroup() %>%
    distinct(term, group) %>%
    sample_n(20)

m1_draws %>%
    ungroup() %>%
    filter(term == "(Intercept)" & str_detect(group, "dep_map_id")) %>%
    mutate(dep_map_id = str_remove(group, "^dep_map_id:")) %>%
    left_join(d %>% distinct(cancer, dep_map_id), by = "dep_map_id") %>%
    arrange(cancer, dep_map_id) %>%
    mutate(dep_map_id = fct_inorder(dep_map_id)) %>%
    ggplot(aes(x = b)) +
    facet_wrap(~ dep_map_id) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_density(aes(color = cancer, fill = cancer), alpha = 0.3, size = 1) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "posterior values",
         y = "probability density",
         title = "Posterior probabilities for varying intercept of cell lines")


m1_draws %>%
    ungroup() %>%
    filter(str_detect(term, "kras_allele")) %>%
    filter(str_detect(group, "hugo_symbol")) %>%
    mutate(kras_allele = str_remove(term, "^kras_allele"),
           info = str_remove(group, "^hugo_symbol:cancer:")) %>%
    separate(info, c("hugo_symbol", "cancer"), sep = ":") %>%
    ggplot(aes(x = kras_allele, y = b)) +
    facet_wrap(cancer ~ hugo_symbol) +
    geom_violin(aes(color = kras_allele, fill = kras_allele)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")

m1_post <- d %>%
    modelr::data_grid(hugo_symbol, cancer, dep_map_id, kras_allele,
                      is_mutated = 0, rna = 0) %>%
    add_predicted_draws(m1)

m1_post %>%
    filter(hugo_symbol %in% c("ATP5PD", "BPHL")) %>%
    group_by(cancer, hugo_symbol, kras_allele) %>%
    median_hdi(.prediction, .width = c(0.5, 0.75, 0.89, 0.95)) %>%
    ggplot(aes(x = kras_allele, y = .prediction, ymin = .lower, ymax = .upper)) +
    geom_pointinterval(aes(color = kras_allele)) +
    facet_wrap(cancer ~ hugo_symbol) +
    scale_color_brewer(palette = "Dark2")


describe_posterior(m1, effects = "all", rope_range = c(-0.1, 0.1)) %>%
    as_tibble() %>%
    arrange(ROPE_Percentage)
