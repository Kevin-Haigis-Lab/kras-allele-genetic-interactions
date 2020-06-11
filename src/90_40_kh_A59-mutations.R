# Analysis of KRAS A59T/E mutants.

GRAPHS_DIR <- "90_40_kh_A59-mutations"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


library(tidybayes)
library(rstanarm)
library(bayestestR)
library(see)



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
#> [1] 38


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
    ) %>%
    knitr::kable(format = "pandoc")
#> is_a59_mutation    frac_in_mapk_kegg_gs   frac_in_mapk_spec_gs
#> ----------------  ---------------------  ---------------------
#> FALSE                         0.9373449              0.6507444
#> TRUE                          0.8823529              0.7058824


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
    ) %>%
    knitr::kable(format = "pandoc")
#> is_a59_mutation    frac_in_mapk_kegg_gs   frac_in_mapk_spec_gs
#> ----------------  ---------------------  ---------------------
#> FALSE                         0.8368486              0.3122415
#> TRUE                          0.8823529              0.6470588


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
    ) %>%
    knitr::kable(format = "pandoc")
#> grp            frac_in_mapk_kegg_gs   frac_in_mapk_spec_gs
#> ------------  ---------------------  ---------------------
#> A59_mutant                0.8823529              0.6470588
#> KRAS_mutant               0.7569480              0.1854472
#> rest                      0.8921946              0.4000700


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
    ) %>%
    knitr::kable(format = "pandoc")
#> grp            frac_in_mapk_kegg_gs   frac_in_mapk_spec_gs
#> ------------  ---------------------  ---------------------
#> A59_mutant                0.8000000              0.6000000
#> KRAS_mutant               0.8847656              0.2285156
#> rest                      0.8767258              0.4201183


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
    ) %>%
    knitr::kable(format = "pandoc")
#> grp            frac_in_mapk_kegg_gs   frac_in_mapk_spec_gs
#> ------------  ---------------------  ---------------------
#> A59_mutant                1.0000000              0.8888889
#> KRAS_mutant               0.9011407              0.5247148
#> rest                      0.9799107              0.8392857


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
#>
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
#>90

# A59 mutants vs. KRAS mutants (KEGG gene set)
d <- a59_mapk_mut_counts %>%
    filter(grp %in% c("A59_mutant", "KRAS_mutant")) %>%
    select(grp, tidyselect::contains("spec")) %>%
    as.data.frame() %>%
    column_to_rownames("grp")
fisher.test(d, alternative = "g")
#>
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
#>


#### ---- Binomial model for impact of A59 mutation ---- ####

is_mapk_mut_ignore_kras <- function(muts, gs) {
    muts_no_kras <- muts[muts != "KRAS"]
    any(muts_no_kras %in% gs)
}

a59_mapk_mut_counts <- a59_cancer_data %>%
    mutate(
        grp = case_when(is_a59_mutation ~ "A59_mutant",
                        ras_allele != "WT" ~ "KRAS_mutant",
                        TRUE ~ "rest"),
        grp = factor(grp, levels = c("rest", "KRAS_mutant", "A59_mutant"))
    ) %>%
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





writeout_model_information <- function(mdl, name) {
    sink(table_path(GRAPHS_DIR, paste0(name, ".txt")))
    for (f in c(summary, describe_posterior, bayestestR::hdi, bf_parameters)) {
        print(f(mdl))
        cat("\n\n")
        cat(str_rep("-", 80), "\n")
        cat("\n\n")
    }
    sink()
}


plot_post_distributions <- function(mdl, name) {
    p <- plot(bayestestR::hdi(mdl, ci = c(0.5, 0.75, 0.89, 0.95)),
              show_intercept = TRUE,
              data = mdl) +
        scale_fill_flat() +
        theme_minimal(base_size = 7, base_family = "Arial") +
        labs(title = glue("89% HDI of {name}"))
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, paste0(name, "_hdi.svg")),
        "medium"
    )
}


#### ---- Fraction of each group that is hypermutated ---- ####

grp_hypermutant_frac <- a59_mapk_mut_counts %>%
    count(grp, is_hypermutant) %>%
    group_by(grp) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
grp_hypermutant_frac %>%
    filter(is_hypermutant) %>%
    knitr::kable(format = "pandoc")
#> grp           is_hypermutant      n        frac
#> ------------  ---------------  ----  ----------
#> rest          TRUE              448   0.1568078
#> KRAS_mutant   TRUE              263   0.1328954
#> A59_mutant    TRUE                9   0.5294118

frac_hypermutant_plot <- grp_hypermutant_frac %>%
    ggplot(aes(x = grp, y = frac)) +
    geom_col(aes(fill = is_hypermutant)) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank()
    ) +
    labs(y = "Fraction of hypermutants",
         title = "Proportion of each group of samples that are hypermutants",
         fill = "is hypermutant")
ggsave_wrapper(
    frac_hypermutant_plot,
    plot_path(GRAPHS_DIR, "frac-hypermutant-plot.svg"),
    "small"
)


#### ---- M1: just `grp` (A59, KRAS, rest) ---- ####

m1_mapk_kegg_stan <- stan_glm(
    is_mapk_kegg_mut ~ 1 + grp,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
writeout_model_information(m1_mapk_kegg_stan, "m1-mapk-kegg-stan")
plot_post_distributions(m1_mapk_kegg_stan, "m1-mapk-kegg-stan")


m1_mapk_spec_stan <- stan_glm(
    is_mapk_spec_mut ~ 1 + grp,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
writeout_model_information(m1_mapk_spec_stan, "m1-mapk-spec-stan")
plot_post_distributions(m1_mapk_spec_stan, "m1-mapk-spec-stan")



#### ---- M2: grp and hypermutant ---- ####

m2_mapk_kegg_stan <- stan_glm(
    is_mapk_kegg_mut ~ 1 + grp + is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
writeout_model_information(m2_mapk_kegg_stan, "m2-mapk-kegg-stan")
plot_post_distributions(m2_mapk_kegg_stan, "m2-mapk-kegg-stan")


m2_mapk_spec_stan <- stan_glm(
    is_mapk_spec_mut ~ 1 + grp + is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
writeout_model_information(m2_mapk_spec_stan, "m2-mapk-spec-stan")
plot_post_distributions(m2_mapk_spec_stan, "m2-mapk-spec-stan")



#### ---- M3: grp and hypermutant and interaction ---- ####


m3_mapk_kegg_stan <- stan_glm(
    is_mapk_kegg_mut ~ 1 + grp * is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
writeout_model_information(m3_mapk_kegg_stan, "m3-mapk-kegg-stan")
plot_post_distributions(m3_mapk_kegg_stan, "m3-mapk-kegg-stan")


m3_mapk_spec_stan <- stan_glm(
    is_mapk_spec_mut ~ 1 + grp * is_hypermutant,
    data = a59_mapk_mut_counts,
    family = binomial(link = "logit"),
    prior_intercept = normal(location = 0, scale = 3),
    prior = normal(location = 0, scale = 3),
    cores = 4,
    seed = 123
)
writeout_model_information(m3_mapk_spec_stan, "m3-mapk-spec-stan")
plot_post_distributions(m3_mapk_spec_stan, "m3-mapk-spec-stan")


#### ---- Compare models ---- ####


m1_mapk_spec_stan$waic <- waic(m1_mapk_spec_stan)
m2_mapk_spec_stan$waic <- waic(m2_mapk_spec_stan)
m3_mapk_spec_stan$waic <- waic(m3_mapk_spec_stan)

m1_mapk_spec_stan$loo <- loo(m1_mapk_spec_stan, cores = 4)
m2_mapk_spec_stan$loo <- loo(m2_mapk_spec_stan, cores = 4)
m3_mapk_spec_stan$loo <- loo(m3_mapk_spec_stan, cores = 4)

model_list <- stanreg_list(m1_mapk_spec_stan, m2_mapk_spec_stan, m3_mapk_spec_stan)
loo_compare(model_list, criterion = "loo")
#> Model comparison based on LOO-CV:
#>                   elpd_diff se_diff
#> m3_mapk_spec_stan    0.0       0.0
#> m2_mapk_spec_stan   -0.4       2.5
#> m1_mapk_spec_stan -310.5      24.6
loo_compare(model_list, criterion = "waic")
#> Model comparison based on WAIC:
#>                   elpd_diff se_diff
#> m3_mapk_spec_stan    0.0       0.0
#> m2_mapk_spec_stan   -0.4       2.5
#> m1_mapk_spec_stan -310.6      24.6

loo_model_weights(model_list)
#> Method: stacking
#> ------
#>                   weight
#> m1_mapk_spec_stan 0.000
#> m2_mapk_spec_stan 0.421
#> m3_mapk_spec_stan 0.579

#### ---- Posterior-Predictive plots ---- ####


post_pred_data <- a59_mapk_mut_counts %>%
    distinct(grp, is_hypermutant) %>%
    arrange(is_hypermutant, grp)
post_pred_data


m3_mapk_spec_post <- post_pred_data %>%
    add_fitted_draws(m3_mapk_spec_stan) %>%
    mean_qi(.width = c(0.89, 0.95))


data_spec_proportions <- a59_mapk_mut_counts %>%
    group_by(grp, is_hypermutant) %>%
    summarise(real_prop = mean(is_mapk_spec_mut)) %>%
    ungroup()

post_pred_plot <- m3_mapk_spec_post %>%
    mutate(condition = paste(grp, is_hypermutant, sep = ", ")) %>%
    ggplot(aes(x = grp)) +
    facet_grid(~ is_hypermutant) +
    geom_linerange(aes(ymin = .lower, ymax = .upper, size = factor(.width))) +
    geom_point(aes(y = .value)) +
    geom_point(aes(y = real_prop),
               data = data_spec_proportions,
               color = "dodgerblue") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_x_discrete(labels = function(x) { str_remove(x, "_mutant") }) +
    scale_size_manual(values = c(1, 0.5)) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank()
    ) +
    labs(title = "Posterior predictions for 'm3_mapk_spec_stan'",
         subtitle = "The separate panels are for non-hypermutants (left) and hypermutants (right).",
         y = "predicted proportion with MAPK mutation",
         size = "PI")
ggsave_wrapper(
    post_pred_plot,
    plot_path(GRAPHS_DIR, "m3-mapk-spec-stan_posterior-pred.svg"),
    "small"
)


sink()