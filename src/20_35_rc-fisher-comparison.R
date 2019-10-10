

test_name_pal <- c(
    RC = "deepskyblue1",
    Fisher = "chocolate1"
)

#### ---- Significant results from each test ---- ####

rc_sig <- rc_test_results %>%
    filter(num_mut_per_cancer >= 3) %>%
    filter(p_val < 0.05 & t_AM >= 3)

fisher_sig <- fisher_comut_df %>%
    filter(p_value_great < 0.05 | p_value_less < 0.05) %>%
    filter(n11 >= 3) %>%
    mutate(test_type = ifelse(
        p_value_great < 0.05, "comutation", "exclusivity"
    ))


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
        kras_allele = str_remove(kras_allele, "KRAS_")) %>%
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

rc_fisher_comparison_barplot_samenum <- bind_rows(rc_sig_samenum, fisher_sig_samenum) %>%
    mutate(
        kras_allele = str_remove(kras_allele, "KRAS_")) %>%
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

