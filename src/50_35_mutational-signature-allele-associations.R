# Test if there is an association between the level of mutational signatures
# and a sample's KRAS allele.

GRAPHS_DIR <- "50_35_mutational-signature-allele-associations"
reset_graph_directory(GRAPHS_DIR)

library(ggasym)
library(ggbeeswarm)

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
    filter(!is_hypermutant) %>%
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


modify_mutsig_names <- function(sig) {
    case_when(
        sig == "N3V2" ~ "N",
        TRUE ~ sig
    )
}

factor_signatures <- function(sig) {
    factor(sig, levels = c(names(mutsig_pal), "N"))
}

plot_allele_sig_association_grids <- function(cancer, data) {
    mod_data <- data %>%
        group_by(signature) %>%
        nest() %>%
        mutate(data = map(data, ~ asymmetrise(.x, group1, group2))) %>%
        unnest(data) %>%
        ungroup() %>%
        mutate(signature = modify_mutsig_names(signature)) %>%
        mutate(signature = factor_signatures(signature),
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
                                guide = guide_colorbar(barwidth = unit(3, "mm"),
                                                       barheight = unit(18, "mm"),
                                                       order = 1)) +
        scale_fill_br_distiller(type = "div",
                                palette = "PiYG",
                                na.value = "grey70",
                                guide = guide_colorbar(barwidth = unit(3, "mm"),
                                                       barheight = unit(18, "mm"),
                                                       order = 2)) +
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

    plt_name <- as.character(glue("{cancer}_allele-sig-associate.svg"))
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, plt_name),
        "large"
    )
    saveFigRds(p, plt_name)
    return(NULL)
}


allele_signature_associations %>%
    group_by(signature, cancer) %>%
    filter(any(p_value < 0.05)) %>%
    select(-data) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(a = walk2(cancer, data, plot_allele_sig_association_grids))


#### ---- Box-plots of specific signatures ---- ####

cancer_sigs_for_boxplot <- tribble(
    ~ cancer, ~ signature,
    "COAD", c("8", "18"),
    "LUAD", c("4", "18"),
    "MM", c("2", "9"),
    "PAAD", c("8", "17")
) %>%
    unnest(signature)


plot_signature_boxplots <- function(cancer, data, min_num_samples = 10) {
    set.seed(0)
    alleles_to_use <- get_alleles_to_test(data, min_num_samples)
    plt_data <- data %>%
        filter(ras_allele %in% alleles_to_use)

    p <- plt_data %>%
        ggplot(aes(ras_allele, contribution)) +
        facet_wrap(~ signature, nrow = 1, scales = "free") +
        geom_boxplot(aes(fill = ras_allele, color = ras_allele),
                     alpha = 0.1, size = 0.5, outlier.shape = NA) +
        geom_quasirandom(aes(color = ras_allele),
                     data = plt_data %>% filter(contribution > 0),
                     size = 0.4, varwidth = TRUE,
                     groupOnX = TRUE, alpha = 0.4) +
        geom_quasirandom(aes(color = ras_allele),
                     data = plt_data %>% filter(contribution == 0),
                     size = 0.4, varwidth = TRUE,
                     groupOnX = TRUE, alpha = 0.4) +
        scale_color_manual(values = short_allele_pal,
                           drop = TRUE, guide = FALSE) +
        scale_fill_manual(values = short_allele_pal,
                          drop = TRUE, guide = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(axis.title.x = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 8, face = "bold")) +
        labs(y = "contribution",
             title = cancer)
    plt_name <- as.character(glue("{cancer}_allele-sig-boxplots.svg"))
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, plt_name),
        "wide"
    )
    saveFigRds(p, plt_name)
    return(NULL)
}

mutsig_noartifact_df %>%
    filter(!is_hypermutant) %>%
    right_join(cancer_sigs_for_boxplot, by = c("cancer", "signature")) %>%
    mutate(ras_allele = str_remove(ras_allele, "KRAS_"),
           ras_allele = factor_alleles(ras_allele),
           signature = factor_signatures(signature)) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(a = walk2(cancer, data, plot_signature_boxplots))
