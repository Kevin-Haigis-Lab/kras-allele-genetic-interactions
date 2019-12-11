
# A palette for some plots in this script.
test_name_pal <- c(
    RC = "deepskyblue1",
    Fisher = "chocolate1"
)

#### ---- Significant results from each test ---- ####

# p-value thresholds for mutual exclusivity and comutation
p_val_cut_mutex <- 0.01
p_val_cut_comut <- 0.01

# Thresholds for the mutation frequency of other gene for mut. ex. and comut.
mutfreq_mutex <- 0.02
mutfreq_comut <- 0.01

# Threshold for the frequency of comutation between the allele and other gene.
comutfreq_comut <- 0.10

# Filter RC-test results using the above thresholds.
rc_sig <- rc_test_results %>%
    filter(
        t_AM >= 10 &
        num_mut_per_cancer / num_samples_per_cancer > !!mutfreq_mutex &
        p_val < !!p_val_cut_mutex
    )

fisher_sig <- fisher_comut_df %>%
    filter(str_replace_us(kras_allele) %in% names(allele_palette)) %>%
    filter(p_value_great < p_val_cut_comut | p_value_less < p_val_cut_mutex) %>%
    filter(
        n11 >= 3 &
        (
            ((n10 + n11) / (n00 + n10 + n01 + n11) > !!mutfreq_comut) |
            (n11 / n01 > !!comutfreq_comut)
        )
    ) %>%
    mutate(test_type = ifelse(
        p_value_great < !!p_val_cut_comut, "comutation", "exclusivity"
    ))



#### ---- Distribution of p-values ---- ####

rc_pval_distribution <- rc_test_results %>%
    filter(num_mut_per_cancer >= 15 & t_AM >= 3) %>%
    ggplot(aes(x = p_val)) +
    facet_grid(rc_test_type ~ cancer, scales = "free") +
    geom_density(aes(color = cancer, fill = cancer), alpha = 0.5) +
    scale_color_manual(values = cancer_palette) +
    scale_fill_manual(values = cancer_palette) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    theme_bw(base_family = "Arial") +
    labs(
        x = "p-value",
        y = "density",
        title = "Distribution of p-values from RC-test"
    )
save_path <- plot_path("20_35_rc-fisher-comparison",
                       "rc_pval_distribution.svg")
ggsave_wrapper(rc_pval_distribution, save_path, "wide")

fish_pval_distribution <- fisher_comut_df %>%
    filter(n10 + n11 > 10) %>%
    ggplot(aes(x = p_value_great)) +
    facet_wrap( ~ cancer, scales = "free") +
    geom_density(aes(color = cancer, fill = cancer), alpha = 0.5) +
    scale_color_manual(values = cancer_palette) +
    scale_fill_manual(values = cancer_palette) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    theme_bw(base_family = "Arial") +
    labs(
        x = "p-value",
        y = "density",
        title = "Distribution of p-values from Fisher's exact test"
    )
save_path <- plot_path("20_35_rc-fisher-comparison",
                       "fish_pval_distribution.svg")
ggsave_wrapper(fish_pval_distribution, save_path, "wide")




#### ---- Comparing gross number of genes ---- ####

rc_summary <- rc_sig %>%
    group_by(cancer, kras_allele, rc_test_type) %>%
    summarise(num_genes_per_type = n_distinct(hugo_symbol)) %>%
    ungroup() %>%
    dplyr::rename(test_type = "rc_test_type") %>%
    mutate(test_name = "RC")

fisher_summary <- fisher_sig %>%
    group_by(cancer, kras_allele, test_type) %>%
    summarise(num_genes_per_type = n_distinct(hugo_symbol)) %>%
    ungroup() %>%
    mutate(test_name = "Fisher")

rc_fisher_comparison_barplot <- bind_rows(rc_summary, fisher_summary) %>%
    mutate(
        kras_allele = str_remove(kras_allele, "KRAS_")
    ) %>%
    ggplot(aes(x = kras_allele, y = num_genes_per_type)) +
    facet_grid(cancer ~ test_type, scales = "free") +
    geom_col(aes(fill = test_name), position = "dodge") +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_fill_manual(values = test_name_pal) +
    theme_classic() +
    theme(
        text = element_text(family = "arial"),
        axis.text.x = element_text(angle = 60, hjust = 1.0),
        legend.position = "bottom"
    )

save_path <- plot_path("20_35_rc-fisher-comparison", "rc_fisher_comparison_barplot.svg")
info(logger, glue("Saving plot {save_path}"))
ggsave_wrapper(rc_fisher_comparison_barplot, save_path, width = 6, height = 8)


# Only looking at the genes with matching number of samples

sample_count_fisher <- fisher_comut_df %>%
    mutate(total_num_fish_samples = n00 + n10 + n01 + n11) %>%
    select(hugo_symbol, cancer, total_num_fish_samples) %>%
    unique()

sample_count_rc <- rc_test_results %>%
    select(hugo_symbol, cancer, num_samples_per_cancer) %>%
    unique()

rc_sig_samenum <- rc_sig %>%
    left_join(sample_count_fisher, by = c("hugo_symbol", "cancer")) %>%
    filter(abs(total_num_fish_samples - num_samples_per_cancer) <= 5) %>%
    group_by(cancer, kras_allele, rc_test_type) %>%
    summarise(num_genes_per_type = n_distinct(hugo_symbol)) %>%
    ungroup() %>%
    dplyr::rename(test_type = "rc_test_type") %>%
    mutate(test_name = "RC")

fisher_sig_samenum <- fisher_sig %>%
    left_join(sample_count_rc, by = c("hugo_symbol", "cancer")) %>%
    filter(abs((n00 + n10 + n01 + n11) - num_samples_per_cancer) <= 5) %>%
    group_by(cancer, kras_allele, test_type) %>%
    summarise(num_genes_per_type = n_distinct(hugo_symbol)) %>%
    ungroup() %>%
    mutate(test_name = "Fisher")

rc_fisher_comparison_barplot_samenum <- bind_rows(
        rc_sig_samenum,
        fisher_sig_samenum
    ) %>%
    mutate(
        kras_allele = str_remove(kras_allele, "KRAS_")
    ) %>%
    ggplot(aes(x = kras_allele, y = num_genes_per_type)) +
    facet_grid(cancer ~ test_type, scales = "free") +
    geom_col(aes(fill = test_name), position = "dodge") +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_fill_manual(values = test_name_pal) +
    theme_classic() +
    theme(
        text = element_text(family = "arial"),
        axis.text.x = element_text(angle = 60, hjust = 1.0),
        legend.position = "bottom"
    )

save_path <- plot_path("20_35_rc-fisher-comparison", "rc_fisher_comparison_barplot_samenum.svg")
info(logger, glue("Saving plot {save_path}"))
ggsave_wrapper(rc_fisher_comparison_barplot_samenum, save_path, width = 6, height = 8)


#### ---- Fisher for Comut. and RC for mut. ex. ---- ####


rc_fisher_comparison_barplot <- bind_rows(
        { rc_summary %>% filter(test_type == "exclusivity") },
        { fisher_summary %>% filter(test_type == "comutation") }
    ) %>%
    mutate(
        kras_allele = str_remove(kras_allele, "KRAS_")
    ) %>%
    ggplot(aes(x = kras_allele, y = num_genes_per_type)) +
    facet_wrap(~ cancer, scales = "free", nrow = 2) +
    geom_col(aes(fill = test_type), position = "dodge") +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_fill_manual(values = comut_mutex_pal) +
    theme_classic() +
    theme(
        text = element_text(family = "arial"),
        axis.text.x = element_text(angle = 40, hjust = 1.0),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank()
    ) +
    labs(
        y = "number of genetic interactions"
    )
save_path <- plot_path("20_35_rc-fisher-comparison", "rc_fisher_comparison_specific.svg")
info(logger, glue("Saving plot {save_path}"))
ggsave_wrapper(rc_fisher_comparison_barplot, save_path, width = 6, height = 8)




#### ---- Combined data frame ---- ####

exclusivity_df <- rc_sig %>%
    filter(rc_test_type == "exclusivity") %>%
    select(
        hugo_symbol, kras_allele, allele, cancer,
        p_val, t_AM,
        num_samples_per_cancer, num_samples_per_cancer_allele,
        num_mut_per_cancer, num_mut_per_cancer_allele
    ) %>%
    mutate(
        test_name = "RC",
        genetic_interaction = "exclusivity"
    )

comutation_df <- fisher_sig %>%
    filter(test_type == "comutation") %>%
    select(
        hugo_symbol, kras_allele, cancer,
        p_value_great, odds_ratio,
        n00, n10, n01, n11,
        allele_freq, gene_freq
    ) %>%
    mutate(
        test_name = "Fisher",
        allele = str_remove(kras_allele, "KRAS_"),
        genetic_interaction = "comutation"
    ) %>%
    dplyr::rename(p_val = p_value_great)

genetic_interaction_df <- bind_rows(exclusivity_df, comutation_df)
ProjectTemplate::cache("genetic_interaction_df")



#### ---- Distribution of the effect size of genetic interactions ---- ####

set.seed(0)
genes_to_label <- c(
    "BRAF", "TP53", "NF1", "NRAS", "PIK3CA",
    unique(cosmic_cgc_df %>% filter(tier == 1) %>% pull(hugo_symbol))
)

rc_mutations_distribition <- genetic_interaction_df %>%
    mutate(
        label = ifelse(
            hugo_symbol %in% !!genes_to_label, hugo_symbol, NA
        ),
        y = ifelse(
            genetic_interaction == "exclusivity",
            t_AM / num_samples_per_cancer,
            n11 / (n00 + n10 + n01 + n11)
        ),
        point_size = ifelse(
            genetic_interaction == "exclusivity",
            num_mut_per_cancer / num_samples_per_cancer,
            gene_freq
        )
    ) %>%
    ggplot(aes(x = -log10(p_val), y = y)) +
    facet_wrap(genetic_interaction ~ cancer, scales = "free", nrow = 2) +
    geom_point(
        aes(alpha = point_size,
            size = point_size,
            color = allele)
    ) +
    ggrepel::geom_text_repel(
        aes(label = label),
        family = "Arial",
        size = 2,
        segment.size = 0.2,
        segment.alpha = 0.5
    ) +
    scale_alpha_continuous(
        range = c(0.1, 0.7),
        guide = guide_legend(label.position = "left", order = 0)
    ) +
    scale_size_continuous(
        range = c(0.1, 3),
        guide = guide_legend(label.position = "left", order = 0)
    ) +
    scale_color_manual(
        values = short_allele_pal,
        guide = guide_legend(label.position = "left", order = 1)
    ) +
    scale_x_continuous(
        expand = expand_scale(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
        expand = expand_scale(mult = c(0.02, 0.02))
    ) +
    theme_bw(
        base_family = "Arial",
        base_size = 8) +
    theme(
        plot.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        legend.key.size = unit(4, "mm")
    ) +
    labs(
        x = "-log10( p-value )",
        y = "number of events per number of samples",
        alpha = "mut. freq",
        size = "mut. freq",
        color = "allele"
    )
save_path <- plot_path("20_35_rc-fisher-comparison",
                       "rc_mutations_distribition.svg")
ggsave_wrapper(rc_mutations_distribition, save_path,
               width = 10, height = 6)
