# Test if there is an association between the level of mutational signatures
# and a sample's KRAS allele.

GRAPHS_DIR <- "50_35_mutational-signature-allele-associations"
reset_graph_directory(GRAPHS_DIR)

library(ggasym)

calc_cohens_d <- function(g1, g2, data) {
    allele_values <- function(a, d) d$contribution[d$ras_allele == a]

    effectsize::cohens_d(x = allele_values(g1, data),
                         y = allele_values(g2, data)) %>%
        as.data.frame() %>%
        pull(Cohens_d) %>%
        unlist()
}


get_alleles_to_test <- function(data, min_num_samples) {
    data %>%
        distinct(tumor_sample_barcode, ras_allele) %>%
        count(ras_allele) %>%
        filter(n >= !!min_num_samples) %>%
        u_pull(ras_allele)
}


pairwise_ttest_msigs <- function(data, min_num_samples = 10) {
    alleles_to_test <- get_alleles_to_test(data, min_num_samples)

    mod_data <- data %>% filter(ras_allele %in% !!alleles_to_test)

    res <- pairwise.wilcox.test(x = mod_data$contribution,
                                g = mod_data$ras_allele,
                                p.adjust.method = "BH",
                                paired = FALSE) %>%
        broom::tidy() %>%
        janitor::clean_names() %>%
        mutate(cohens_d = map2_dbl(group1, group2, calc_cohens_d, data = data))

    return(res)
}


allele_signature_associations <- mutsig_noartifact_df %>%
    group_by(cancer, signature) %>%
    nest() %>%
    mutate(t_test_res = map(data, pairwise_ttest_msigs)) %>%
    unnest(t_test_res) %>%
    ungroup()


assign_stars <- function(pval) {
    case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        pval < 0.05 ~ "*",
        TRUE ~ NA_character_
    )
}


plot_allele_sig_association_grids <- function(cancer, data) {
    mod_data <- data %>%
        group_by(signature) %>%
        nest() %>%
        mutate(data = map(data, ~ asymmetrise(.x, group1, group2))) %>%
        unnest(data) %>%
        ungroup() %>%
        mutate(signature = factor(signature, levels = names(mutsig_pal)),
               group1 = str_remove(group1, "KRAS_"),
               group2 = str_remove(group2, "KRAS_"),
               group1 = factor_alleles(group1, reverse = FALSE),
               group2 = factor_alleles(group2, reverse = FALSE))

    p <- mod_data %>%
        mutate(p_value = ifelse(is.na(p_value), 1, p_value),
               cohens_d = ifelse(is.na(cohens_d), 0, cohens_d),
               stars = assign_stars(p_value)) %>%
        ggplot(aes(x = group1, y = group2)) +
        facet_wrap(~ signature) +
        geom_asymmat(aes(fill_tl = p_value, fill_br = cohens_d),
                     fill_diag = "grey80", color = "white") +
        geom_text(aes(label = stars), size = 2) +
        scale_fill_tl_distiller(type = "seq",
                                palette = "OrRd",
                                na.value = "grey70",
                                guide = guide_colorbar(order = 1)) +
        scale_fill_br_distiller(type = "div",
                                palette = "PiYG",
                                na.value = "grey70",
                                guide = guide_colorbar(order = 2)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(axis.title = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 8, face = "bold"),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(fill_br = "Cohen's d",
             fill_tl = "adj. p-value",
             title = cancer)
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, glue("{cancer}_allele-sig-associate.svg")),
        "large"
    )
    return(NULL)
}


allele_signature_associations %>%
    group_by(signature, cancer) %>%
    filter(any(p_value < 0.05)) %>%
    select(-data) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(a = walk2(cancer, data, plot_allele_sig_association_grids))
