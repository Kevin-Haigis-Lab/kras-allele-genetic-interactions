# A second method for obtaining KRAS allele-specific synthetic lethal targets
# for Shikha.


library(tidybayes)
library(bayesplot)
library(easystats)
library(rstanarm)



GRAPHS_DIR <- "90_31_synlet-for-shikha_m2"
reset_graph_directory(GRAPHS_DIR)

TABLES_DIR <- GRAPHS_DIR
reset_table_directory(TABLES_DIR)

set.seed(0)

#### ---- Data preparation ---- ####

# Isolate data that is pertinent to this analysis.
modeling_data <- depmap_modelling_df %>%
  filter_depmap_by_allele_count() %>%
  group_by(cancer) %>%
  filter(n_distinct(kras_allele) >= 3) %>%
  filter(!is_deleted) %>%
  add_count(cancer, hugo_symbol, kras_allele) %>%
  group_by(cancer, hugo_symbol) %>%
  filter(all(n >= 3)) %>%
  select(-n) %>%
  ungroup()

# Only COAD and scale RNA expression.
modeling_data <- modeling_data %>%
  filter(cancer == "COAD") %>%
  group_by(hugo_symbol) %>%
  mutate(rna_expression_std = scale_numeric(rna_expression,
    na.rm = TRUE
  )) %>%
  ungroup()


#### ---- Data for experimentation ---- ####

set.seed(0)
eg_genes <- sample(unique(modeling_data$hugo_symbol), 10)
eg_genes <- c(
  eg_genes,
  "STARD9",
  "PAX6",
  "KNTC1",
  "FAF2",
  "THSD7A",
  "SERPING1",
  "LIN7C",
  "DCHS1",
  "ORC4",
  "MTRF1L"
)
n_distinct(eg_genes)
d <- modeling_data %>%
  filter(hugo_symbol %in% eg_genes)

d %>%
  write_tsv(table_path(TABLES_DIR, "sample-modeling-data.tsv"))

## NOTE: Model experimentation in "90_33_synlet-for-shikha_m2-expt.R"


#### ---- Modeling ---- ####

all_kras_alleles <- sort(unique(modeling_data$kras_allele))
kras_allele_order <- c("WT", all_kras_alleles[all_kras_alleles != "WT"])

modeling_data <- modeling_data %>%
  select(hugo_symbol:dep_map_id, gene_effect:kras_allele,
    is_mutated,
    rna = rna_expression_std
  ) %>%
  mutate(
    kras_allele = factor(kras_allele, levels = kras_allele_order),
    is_mutated = as.numeric(is_mutated)
  )



ProjectTemplate::cache("coad_kras_sl_m1", depends = "modeling_data", {
  coad_kras_sl_m1 <- stan_glmer(
    gene_effect ~ 1 + (1 + rna + kras_allele + is_mutated | hugo_symbol),
    data = modeling_data,
    family = gaussian(link = "identity"),
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    refresh = 100,
    cores = 4
  )

  coad_kras_sl_m1$loo <- loo(coad_kras_sl_m1, cores = 4)

  return(coad_kras_sl_m1)
})


ProjectTemplate::cache("coad_kras_sl_m2", depends = "modeling_data", {
  coad_kras_sl_m2 <- stan_glmer(
    gene_effect ~ 1 + (1 + rna + kras_allele | hugo_symbol),
    data = modeling_data,
    family = gaussian(link = "identity"),
    prior = student_t(df = 2, location = 0, scale = 1),
    prior_intercept = normal(location = 0, scale = 0.5),
    prior_aux = cauchy(),
    prior_covariance = decov(),
    chains = 4,
    adapt_delta = 0.99,
    iter = 4000,
    refresh = 100,
    cores = 4
  )

  coad_kras_sl_m2$loo <- loo(coad_kras_sl_m2, cores = 4)

  return(coad_kras_sl_m2)
})
