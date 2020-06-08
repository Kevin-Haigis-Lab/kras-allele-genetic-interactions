# Analysis of KRAS A59T/E mutants.

GRAPHS_DIR <- "90_40_kh_A59-mutations"
reset_graph_directory(GRAPHS_DIR)

a59_muts <- cancer_full_coding_muts_df %>%
    filter(cancer == "COAD") %>%
    filter(hugo_symbol == "KRAS") %>%
    filter(amino_acid_change %in% c("A59E", "A59T")) %>%
    mutate(seq_type = case_when(
        target == "exome" ~ "exome",
        target == "genome" ~ "genome",
        TRUE ~ "targeted"
    )) %>%
    select(tumor_sample_barcode,
           a59_mutation = amino_acid_change, seq_type)

a59_muts %>%
    count(a59_mutation, seq_type) %>%
    arrange(seq_type, -n)

mapk_kegg_gs <- kegg_geneset_df %>%
    filter(gene_set == "Mapk signaling pathway") %>%
    u_pull(hugo_symbol)

mapk_spec_gs <- c("KRAS", "NF1", "BRAF", "RAF1", "ARAF", "NRAS", "HRAS",
                     "MAP2K1", "EGFR", "ERBB2", "ERBB3", "ERBB4", "RASAL1",
                     "RASAL2")

# A specific data frame for this analysis (only in COAD).
a59_cancer_data <- cancer_full_coding_muts_df %>%
    filter(cancer == "COAD") %>%
    left_join(a59_muts %>% select(tumor_sample_barcode, a59_mutation),
              by = "tumor_sample_barcode") %>%
    mutate(a59_mutation = ifelse(is.na(a59_mutation), "other", a59_mutation),
           is_a59_mutation = a59_mutation != "other",
           in_mapk_kegg_gs = hugo_symbol %in% !!mapk_kegg_gs,
           in_mapk_spec_gs = hugo_symbol %in% !!mapk_spec_gs)


# Any tumor samples lost by removing KRAS mutations?
# These samples would only have a single mutation and it would be a KRAS mut.
num_samples <- n_distinct(a59_cancer_data$tumor_sample_barcode)
num_samples_woKRAS <- a59_cancer_data %>%
    filter(hugo_symbol != "KRAS") %>%
   pull(tumor_sample_barcode) %>%
    n_distinct()
# Number of samples with only KRAS mutations.
num_samples - num_samples_woKRAS


#### ---- A59 muts vs. rest freq. of MAPK muts without KRAS ---- ####

num_samples_df <- a59_cancer_data %>%
    distinct(tumor_sample_barcode, is_a59_mutation) %>%
    count(is_a59_mutation)

a59_cancer_data %>%
    filter(
        !(hugo_symbol == "KRAS" & amino_acid_change %in% c("A59E", "A59T"))
    ) %>%
    group_by(tumor_sample_barcode, is_a59_mutation) %>%
    summarise(in_mapk_kegg_gs = any(in_mapk_kegg_gs),
              in_mapk_spec_gs = any(in_mapk_spec_gs)) %>%
    group_by(is_a59_mutation) %>%
    summarise(num_in_mapk_kegg_gs = sum(in_mapk_kegg_gs),
              num_in_mapk_spec_gs = sum(in_mapk_spec_gs)) %>%
    ungroup() %>%
    left_join(num_samples_df, by = "is_a59_mutation") %>%
    transmute(
        is_a59_mutation,
        frac_in_mapk_kegg_gs = num_in_mapk_kegg_gs / n,
        frac_in_mapk_spec_gs = num_in_mapk_spec_gs / n
    )


#### ---- A59 muts vs. rest freq. of MAPK muts without KRAS ---- ####

num_samples_df <- a59_cancer_data %>%
    distinct(tumor_sample_barcode, is_a59_mutation) %>%
    count(is_a59_mutation)

a59_cancer_data %>%
    filter(hugo_symbol != "KRAS") %>%
    group_by(tumor_sample_barcode, is_a59_mutation) %>%
    summarise(in_mapk_kegg_gs = any(in_mapk_kegg_gs),
              in_mapk_spec_gs = any(in_mapk_spec_gs)) %>%
    group_by(is_a59_mutation) %>%
    summarise(num_in_mapk_kegg_gs = sum(in_mapk_kegg_gs),
              num_in_mapk_spec_gs = sum(in_mapk_spec_gs)) %>%
    ungroup() %>%
    left_join(num_samples_df, by = "is_a59_mutation") %>%
    transmute(
        is_a59_mutation,
        frac_in_mapk_kegg_gs = num_in_mapk_kegg_gs / n,
        frac_in_mapk_spec_gs = num_in_mapk_spec_gs / n
    )


#### ---- A59 muts vs. KRAS muts. vs. rest freq. of MAPK muts ---- ####

num_samples_df <- a59_cancer_data %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    distinct(tumor_sample_barcode, grp) %>%
    count(grp)

a59_cancer_data %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    filter(hugo_symbol != "KRAS") %>%
    group_by(tumor_sample_barcode, grp) %>%
    summarise(in_mapk_kegg_gs = any(in_mapk_kegg_gs),
              in_mapk_spec_gs = any(in_mapk_spec_gs)) %>%
    group_by(grp) %>%
    summarise(num_in_mapk_kegg_gs = sum(in_mapk_kegg_gs),
              num_in_mapk_spec_gs = sum(in_mapk_spec_gs)) %>%
    ungroup() %>%
    left_join(num_samples_df, by = "grp") %>%
    transmute(
        grp,
        frac_in_mapk_kegg_gs = num_in_mapk_kegg_gs / n,
        frac_in_mapk_spec_gs = num_in_mapk_spec_gs / n
    )


#### ---- A59 muts vs. KRAS muts. vs. rest freq. of MAPK muts (WG/ES) ---- ####

num_samples_df <- a59_cancer_data %>%
    filter(target %in% c("genome", "exome")) %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    distinct(tumor_sample_barcode, grp) %>%
    count(grp)

a59_cancer_data %>%
    filter(target %in% c("genome", "exome")) %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    filter(hugo_symbol != "KRAS") %>%
    group_by(tumor_sample_barcode, grp) %>%
    summarise(in_mapk_kegg_gs = any(in_mapk_kegg_gs),
              in_mapk_spec_gs = any(in_mapk_spec_gs)) %>%
    group_by(grp) %>%
    summarise(num_in_mapk_kegg_gs = sum(in_mapk_kegg_gs),
              num_in_mapk_spec_gs = sum(in_mapk_spec_gs)) %>%
    ungroup() %>%
    left_join(num_samples_df, by = "grp") %>%
    transmute(
        grp,
        frac_in_mapk_kegg_gs = num_in_mapk_kegg_gs / n,
        frac_in_mapk_spec_gs = num_in_mapk_spec_gs / n
    )


#### ---- A59 muts vs. KRAS muts. vs. rest freq. of MAPK muts (hypermuts) ---- ####

num_samples_df <- a59_cancer_data %>%
    filter(is_hypermutant) %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    distinct(tumor_sample_barcode, grp) %>%
    count(grp)

a59_cancer_data %>%
    filter(is_hypermutant) %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    filter(hugo_symbol != "KRAS") %>%
    group_by(tumor_sample_barcode, grp) %>%
    summarise(in_mapk_kegg_gs = any(in_mapk_kegg_gs),
              in_mapk_spec_gs = any(in_mapk_spec_gs)) %>%
    group_by(grp) %>%
    summarise(num_in_mapk_kegg_gs = sum(in_mapk_kegg_gs),
              num_in_mapk_spec_gs = sum(in_mapk_spec_gs)) %>%
    ungroup() %>%
    left_join(num_samples_df, by = "grp") %>%
    transmute(
        grp,
        frac_in_mapk_kegg_gs = num_in_mapk_kegg_gs / n,
        frac_in_mapk_spec_gs = num_in_mapk_spec_gs / n
    )


#### ---- Fisher of MAPK mutation between A59 muts and KRAS muts ---- ####

num_samples_df <- a59_cancer_data %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    distinct(tumor_sample_barcode, grp) %>%
    count(grp)

a59_mapk_mut_counts <- a59_cancer_data %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    filter(hugo_symbol != "KRAS") %>%
    group_by(tumor_sample_barcode, grp) %>%
    summarise(in_mapk_kegg_gs = any(in_mapk_kegg_gs),
              in_mapk_spec_gs = any(in_mapk_spec_gs)) %>%
    group_by(grp) %>%
    summarise(num_in_mapk_kegg_gs = sum(in_mapk_kegg_gs),
              num_in_mapk_spec_gs = sum(in_mapk_spec_gs)) %>%
    ungroup() %>%
    left_join(num_samples_df, by = "grp") %>%
    mutate(num_not_mapk_kegg_gs = n - num_in_mapk_kegg_gs,
           num_not_mapk_spec_gs = n - num_in_mapk_spec_gs) %>%
    select(grp, tidyselect::starts_with("num"), total = n)

# A59 mutants vs. KRAS mutants (KEGG gene set)
d <- a59_mapk_mut_counts %>%
    filter(grp %in% c("A59_mutant", "KRAS_mutant")) %>%
    select(grp, tidyselect::contains("kegg")) %>%
    as.data.frame() %>%
    column_to_rownames("grp")
fisher.test(d, alternative = "g")
#>     Fisher's Exact Test for Count Data
#>
#> data:  d
#> p-value = 0.1815
#> alternative hypothesis: true odds ratio is greater than 1
#> 95 percent confidence interval:
#>  0.6598959       Inf
#> sample estimates:
#> odds ratio
#>   2.407507

# A59 mutants vs. KRAS mutants (KEGG gene set)
d <- a59_mapk_mut_counts %>%
    filter(grp %in% c("A59_mutant", "KRAS_mutant")) %>%
    select(grp, tidyselect::contains("spec")) %>%
    as.data.frame() %>%
    column_to_rownames("grp")
fisher.test(d, alternative = "g")
#>     Fisher's Exact Test for Count Data
#>
#> data:  d
#> p-value = 4.063e-05
#> alternative hypothesis: true odds ratio is greater than 1
#> 95 percent confidence interval:
#>  3.150826      Inf
#> sample estimates:
#> odds ratio
#>   8.041764


#### ---- Binomial model for impact of A59 mutation ---- ####

is_mapk_mut_ignore_kras <- function(muts, gs) {
    muts_no_kras <- muts[muts != "KRAS"]
    any(muts_no_kras %in% gs)
}

a59_mapk_mut_counts <- a59_cancer_data %>%
    mutate(grp = case_when(is_a59_mutation ~ "A59_mutant",
                           ras_allele != "WT" ~ "KRAS_mutant",
                           TRUE ~ "rest")) %>%
    group_by(tumor_sample_barcode, is_hypermutant, grp) %>%
    summarise(
        is_mapk_kegg_mut = is_mapk_mut_ignore_kras(hugo_symbol,
                                                   gs = mapk_kegg_gs),
        is_mapk_kegg_mut = as.numeric(is_mapk_kegg_mut),
        is_mapk_spec_mut = is_mapk_mut_ignore_kras(hugo_symbol,
                                                   gs = mapk_spec_gs),
        is_mapk_spec_mut = as.numeric(is_mapk_spec_mut)
    ) %>%
    ungroup()

a59_mapk_mut_counts %>%
    group_by(grp) %>%
    summarise(frac_mapk_kegg_mut = sum(is_mapk_kegg_mut) / n(),
            frac_mapk_spec_mut = sum(is_mapk_spec_mut) / n())



m_mapk_kegg <- glm(is_mapk_kegg_mut ~ grp,
                   data = a59_mapk_mut_counts,
                   family = "binomial")
summary(m_mapk_kegg)
exp(coef(m_mapk_kegg))

m_mapk_spec <- glm(is_mapk_spec_mut ~ grp,
                   data = a59_mapk_mut_counts,
                   family = "binomial")
summary(m_mapk_spec)
exp(coef(m_mapk_spec))


library(rstanarm)
library(bayestestR)
library(see)

m1_mapk_kegg_stan <- stan_glm(
    is_mapk_kegg_mut ~ 1 + grp,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
summary(m1_mapk_kegg_stan)
describe_posterior(m1_mapk_kegg_stan)
hdi(m1_mapk_kegg_stan)

m1_mapk_kegg_stan_hpi_plot <- plot(hdi(m1_mapk_kegg_stan),
                                   show_intercept = TRUE) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    labs(title = "(model 1) KRAS mutation group")
ggsave_wrapper(
    m1_mapk_kegg_stan_hpi_plot,
    plot_path(GRAPHS_DIR, "m1_mapk_kegg_stan_hpi_plot.svg"),
    "medium"
)



m2_mapk_kegg_stan <- stan_glm(
    is_mapk_kegg_mut ~ 1 + grp + is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
summary(m2_mapk_kegg_stan)
describe_posterior(m2_mapk_kegg_stan)
hdi(m2_mapk_kegg_stan)

m2_mapk_kegg_stan_hpi_plot <- plot(hdi(m2_mapk_kegg_stan),
                                   show_intercept = TRUE) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    labs(title = "(model 2) KRAS mutation group and hypermutant")
ggsave_wrapper(
    m2_mapk_kegg_stan_hpi_plot,
    plot_path(GRAPHS_DIR, "m2_mapk_kegg_stan_hpi_plot.svg"),
    "medium"
)

post_pred_data <- a59_mapk_mut_counts %>%
    distinct(grp, is_hypermutant) %>%
    arrange(is_hypermutant, grp)
m2_post <- posterior_predict(m2_mapk_kegg_stan, newdata = post_pred_data)

post_pred_data %>%
    mutate(avg_pred = apply(m2_post, 2, mean))




m2_mapk_spec_stan <- stan_glm(
    is_mapk_spec_mut ~ 1 + grp + is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
summary(m2_mapk_spec_stan)
describe_posterior(m2_mapk_spec_stan)
hdi(m2_mapk_spec_stan)

m2_mapk_spec_stan_hpi_plot <- plot(hdi(m2_mapk_spec_stan),
                                   show_intercept = TRUE) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    labs(title = "(model 2) KRAS mutation group and hypermutant")
ggsave_wrapper(
    m2_mapk_spec_stan_hpi_plot,
    plot_path(GRAPHS_DIR, "m2_mapk_spec_stan_hpi_plot.svg"),
    "medium"
)





m3_mapk_kegg_stan <- stan_glm(
    is_mapk_kegg_mut ~ 1 + grp*is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
summary(m3_mapk_kegg_stan)
describe_posterior(m3_mapk_kegg_stan)
hdi(m3_mapk_kegg_stan)

m3_mapk_kegg_stan_hpi_plot <- plot(hdi(m3_mapk_kegg_stan),
                                   show_intercept = TRUE) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    labs(title = "(model 3) KRAS mutation group and hypermutant and interaction")
ggsave_wrapper(
    m3_mapk_kegg_stan_hpi_plot,
    plot_path(GRAPHS_DIR, "m3_mapk_kegg_stan_hpi_plot.svg"),
    "medium"
)



m3_mapk_spec_stan <- stan_glm(
    is_mapk_spec_mut ~ 1 + grp*is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
summary(m3_mapk_spec_stan)
describe_posterior(m3_mapk_spec_stan)
hdi(m3_mapk_spec_stan)

m3_mapk_spec_stan_hpi_plot <- plot(hdi(m3_mapk_spec_stan),
                                   show_intercept = TRUE) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    labs(title = "(model 3) KRAS mutation group and hypermutant and interaction")
ggsave_wrapper(
    m3_mapk_spec_stan_hpi_plot,
    plot_path(GRAPHS_DIR, "m3_mapk_spec_stan_hpi_plot.svg"),
    "medium"
)