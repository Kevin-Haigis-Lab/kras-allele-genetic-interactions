

library(mustashe)
library(jhcutils)
library(glue)
library(ggtext)
library(patchwork)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern ="R$")) {
    if (str_detect(f, "global|enrich")) { next }
    source(f)
}

options(dplyr.summarise.inform = FALSE)

real_kras_mutations <- readRDS("~/Downloads/real_kras_mutations.rds")
kras_allele_predictions <- readRDS("~/Downloads/kras_allele_predictions.rds")
all_kras_allele_predictions <- readRDS("~/Downloads/all_kras_allele_predictions.rds")


theme_set(theme_bw(base_size = 11, base_family = "Arial"))

ranked_allele_predictions <- kras_allele_predictions %>%
    inner_join(real_kras_mutations %>% rename(real_kras_allele = kras_allele),
               by = c("cancer", "tumor_sample_barcode")) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    arrange(-allele_prob) %>%
    mutate(allele_idx = row_number()) %>%
    ungroup() %>%
    arrange(cancer, tumor_sample_barcode, allele_idx)



ranked_allele_predictions %>%
    filter(real_kras_allele == kras_allele) %>%
    count(cancer, real_kras_allele, allele_idx) %>%
    mutate(top_2 = allele_idx < 3) %>%
    ggplot(aes(real_kras_allele, n)) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_col(aes(fill = fct_rev(factor(allele_idx)))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_fill_brewer(type = "div", palette = "RdYlBu") +
    theme(strip.background = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(x = NULL,
         y = "number of tumor samples",
         fill = "rank",
         title = "The rank of the observed allele when ordered by probability using mutational signatures")


ranked_allele_predictions %>%
    filter(real_kras_allele == "G12C" & cancer == "LUAD") %>%
    filter(kras_allele != "G13C") %>%
    ggplot(aes(x = allele_prob)) +
    geom_density(aes(color = kras_allele, fill = kras_allele), alpha = 0.1)


ranked_probability_plot <- function(df, allele, cancer, ignore_alleles = c()) {
    title <- glue("The distribution of probabilities in KRAS {allele} mutant {cancer}")
    df %>%
        filter(real_kras_allele == !!allele & cancer %in% !!cancer) %>%
        filter(!(kras_allele %in% !!ignore_alleles)) %>%
        mutate(kras_allele = fct_reorder(kras_allele, allele_prob)) %>%
        group_by(tumor_sample_barcode) %>%
        mutate(hit = any(allele_idx == 1 & kras_allele == real_kras_allele)) %>%
        ungroup() %>%
        mutate(hit = factor(hit, levels = c("TRUE", "FALSE"))) %>%
        ggplot(aes(x = kras_allele, y = allele_prob)) +
        geom_line(aes(group = tumor_sample_barcode, color = hit, alpha = hit)) +
        geom_point(aes(color = hit), alpha = 0.4, size = 3) +
        scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                           limits = c(0, NA)) +
        scale_color_manual(values = c("red", "black")) +
        scale_alpha_manual(values = c(0.4, 0.3)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(x = "possible KRAS allele",
             y = "probability of allele from mutational signatures",
             title = title)
}


ranked_allele_predictions %>%
    ranked_probability_plot("G12C", "LUAD", "G13C")

ranked_allele_predictions %>%
    ranked_probability_plot("G12R", "PAAD")

ranked_allele_predictions %>%
    ranked_probability_plot("G13D", "COAD")

ranked_allele_predictions %>%
    ranked_probability_plot("G12D", "COAD")

ranked_allele_predictions %>%
    ranked_probability_plot("A146T", "COAD")
