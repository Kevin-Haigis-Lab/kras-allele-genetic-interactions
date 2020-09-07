# Test if there is an association between the level of mutational signatures
# and a sample's KRAS allele.

GRAPHS_DIR <- "50_35_mutational-signature-allele-associations"
reset_graph_directory(GRAPHS_DIR)

library(ggasym)
library(ggbeeswarm)
library(ggridges)

calc_cohens_d <- function(g1, g2, data) {
    allele_values <- function(a, d) d$contribution[d$ras_allele == a]

    effectsize::cohens_d(x = allele_values(g1, data),
                         y = allele_values(g2, data)) %>%
        as.data.frame() %>%
        pull(Cohens_d) %>%
        unlist()
}


calc_diff_medians <- function(g1, g2, data) {
    allele_values <- function(a, d) d$contribution[d$ras_allele == a]
    median(allele_values(g1, data)) - median(allele_values(g2, data))
}


get_alleles_to_test <- function(data, min_num_samples) {
    data %>%
        distinct(tumor_sample_barcode, ras_allele) %>%
        count(ras_allele) %>%
        filter(n >= !!min_num_samples) %>%
        u_pull(ras_allele)
}


pairwise_ttest_msigs <- function(data) {
    res <- pairwise.wilcox.test(x = data$contribution,
                                g = data$ras_allele,
                                p.adjust.method = "BH",
                                paired = FALSE) %>%
        broom::tidy() %>%
        janitor::clean_names() %>%
        mutate(cohens_d = map2_dbl(group1, group2,
                                   calc_cohens_d, data = data),
               diff_med = map2_dbl(group1, group2,
                                   calc_diff_medians, data = data))

    return(res)
}


pairwise_ttest_cancer_msigs <- function(data, min_num_samples = 10) {
    alleles_to_test <- get_alleles_to_test(data, min_num_samples)
    mod_data <- data %>% filter(ras_allele %in% !!alleles_to_test)
    pairwise_ttest_msigs(mod_data)
}



allele_signature_associations <- mutsig_noartifact_df %>%
    filter(!is_hypermutant) %>%
    group_by(cancer, signature) %>%
    nest() %>%
    mutate(t_test_res = map(data, pairwise_ttest_cancer_msigs)) %>%
    unnest(t_test_res) %>%
    ungroup()


modify_mutsig_names <- function(sig) {
    case_when(
        sig == "N3V2" ~ "N",
        TRUE ~ sig
    )
}


strip_ras <- function(df, col = ras_allele) {
    df %>%
        mutate({{col}} := str_remove({{col}}, "KRAS_"))
}


factor_signatures <- function(sig) {
    factor(sig, levels = c(names(mutsig_pal), "N"))
}


plot_sig_associations_ggasym <- function(data,
                                         fill_br, fill_br_lbl,
                                         title = NULL) {
    data %>%
        mutate(p_value = ifelse(is.na(p_value), 1, p_value),
               fill_br = ifelse(is.na({{ fill_br }}), 0, {{ fill_br }}),
               stars = assign_stars(p_value)) %>%
        ggplot(aes(x = group1, y = group2)) +
        geom_asymmat(aes(fill_tl = p_value, fill_br = fill_br),
                     fill_diag = "grey80", color = "white") +
        geom_text(aes(label = stars), size = 2) +
        scale_fill_tl_distiller(
            type = "seq",
            palette = "OrRd",
            na.value = "grey70",
            guide = guide_colorbar(barwidth = unit(3, "mm"),
                                   barheight = unit(18, "mm"),
                                   order = 1)) +
        scale_fill_br_distiller(
            type = "div",
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
        labs(fill_br = fill_br_lbl,
             fill_tl = "adj. p-value",
             title = title)
}


plot_allele_sig_association_grids <- function(cancer,
                                              data,
                                              fill_br,
                                              fill_br_lbl = NULL,
                                              plt_name_glue) {
    mod_data <- data %>%
        group_by(signature) %>%
        nest() %>%
        mutate(data = map(data, ~ asymmetrise(.x, group1, group2))) %>%
        unnest(data) %>%
        ungroup() %>%
        strip_ras(group1) %>%
        strip_ras(group2) %>%
        mutate(signature = modify_mutsig_names(signature),
               signature = factor_signatures(signature),
               group1 = factor_alleles(group1, reverse = FALSE),
               group2 = factor_alleles(group2, reverse = FALSE))


    p <- mod_data %>%
        plot_sig_associations_ggasym(fill_br = {{ fill_br }},
                                     fill_br_lbl = fill_br_lbl,
                                     title = cancer) +
        facet_wrap(~ signature)

    plt_name <- as.character(glue(plt_name_glue))
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
    mutate(
        a = walk2(cancer, data, plot_allele_sig_association_grids,
                  fill_br = cohens_d,
                  fill_br_lbl = "Cohen's d",
                  plt_name_glue = "{cancer}_cohensd_allele-sig-associate.svg"),
        a = walk2(cancer, data, plot_allele_sig_association_grids,
                  fill_br = diff_med,
                  fill_br_lbl = "diff. in median",
                  plt_name_glue = "{cancer}_median_allele-sig-associate.svg")
    )


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
    strip_ras(ras_allele) %>%
    mutate(ras_allele = factor_alleles(ras_allele),
           signature = factor_signatures(signature)) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(a = walk2(cancer, data, plot_signature_boxplots))


################################################################################
################################################################################
##> Spare plots for PJP


#### ---- Mut. Sig. 4 in LUAD as an example ---- ####

mutsig4_luad_df <- mutsig_noartifact_df %>%
    filter(cancer == "LUAD" & signature == "4")


mutsig4_luad_tiles <- mutsig4_luad_df %>%
    group_by(ras_allele) %>%
    filter(n_distinct(tumor_sample_barcode) > 10) %>%
    ungroup() %>%
    arrange(ras_allele, contribution) %>%
    mutate(idx = row_number()) %>%
    strip_ras() %>%
    make_rows_cols_idx(idx, 21) %>%
    ggplot(aes(x = row_idx, y = col_idx)) +
    geom_tile(aes(fill = ras_allele, alpha = contribution),
              color = "grey20") +
    coord_fixed() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = short_allele_pal, drop = TRUE) +
    theme(legend.title = element_markdown(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(3, "mm")) +
    labs(x = NULL,
         y = NULL,
         fill = "*KRAS* allele",
         contribution = "signature level",
         title = "Mutational signature 4 levels per sample in LUAD",
         subtitle = "Each cell represents an individual tumor sample.")
ggsave_wrapper(mutsig4_luad_tiles,
               plot_path(GRAPHS_DIR, "mutsig4-luad_tiles.svg"),
               "small")


mutsig4_luad_bar <- mutsig4_luad_df %>%
    group_by(ras_allele) %>%
    filter(n_distinct(tumor_sample_barcode) > 10) %>%
    ungroup() %>%
    arrange(contribution) %>%
    strip_ras() %>%
    mutate(tumor_sample_barcode = fct_inorder(tumor_sample_barcode),
           ras_allele = fct_reorder(ras_allele, contribution)) %>%
    ggplot(aes(x = ras_allele)) +
    geom_bar(aes(alpha = contribution, group = contribution),
             width = 1.0, fill = "black", position = "fill") +
    scale_alpha_continuous(range = c(0.05, 1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.title.x = element_markdown(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    labs(x = "*KRAS* allele of tumor",
         y = "individual tumor samples",
         alpha = "mutational signature 4 level",
         title = "Mutational signature 4 level in LUAD tumor samples")
ggsave_wrapper(mutsig4_luad_bar,
               plot_path(GRAPHS_DIR, "mutsig4-luad_bar.svg"),
               "small")


mutsig_4_ggridge <- mutsig4_luad_df %>%
    group_by(ras_allele) %>%
    filter(n_distinct(tumor_sample_barcode) > 10) %>%
    ungroup() %>%
    mutate(ras_allele = str_remove(ras_allele, "KRAS_"),
           ras_allele = fct_reorder(ras_allele, contribution)) %>%
    ggplot(aes(x = contribution, y = ras_allele)) +
    geom_density_ridges(aes(color = ras_allele, fill = ras_allele),
                        size = 1, alpha = 0.5) +
    scale_color_manual(values = short_allele_pal, drop = TRUE) +
    scale_fill_manual(values = short_allele_pal, drop = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(panel.background = element_blank(),
          legend.position = "none",
          plot.title = element_markdown(hjust = 0.5),
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown()) +
    labs(x = "mutational signature 4 levels in tumor samples",
         y = "distribution density",
         title = "Distribution of signature 4 levels in LUAD tumor samples when<br>separated by their
         *KRAS* mutation")
ggsave_wrapper(mutsig_4_ggridge,
               plot_path(GRAPHS_DIR, "mutsig-4-ggridge.svg"),
               "small")


#### ---- All mut. sigs. ---- ####

plot_memosorted_tumor_samples <- function(cancer, data) {
    p <- data %>%
        group_by(ras_allele) %>%
        filter(n_distinct(tumor_sample_barcode) > 10) %>%
        mutate(description = factor(description, levels = names(mutsig_descrpt_pal)),
               score = as.numeric(fct_rev(description)) ** 2 * contribution,
               tsb_idx = as.numeric(fct_reorder(tumor_sample_barcode,
                                                score,
                                                .fun = sum)),
               ras_allele = factor_alleles(ras_allele)) %>%
        ggplot(aes(x = tsb_idx, y = contribution)) +
        facet_grid(~ ras_allele, scales = "free_x", space = "free") +
        geom_col(aes(fill = description), width = 1.0) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_fill_manual(values = mutsig_descrpt_pal,
                          drop = TRUE,
                          guide = guide_legend(nrow = 1,
                                               label.position = "top")) +
        theme(axis.text.x = element_blank(),
              strip.background = element_blank(),
              strip.text = element_markdown(angle = 90),
              legend.position = "bottom",
              legend.key.size = unit(3, "mm"),
              plot.title = element_markdown(hjust = 0.5)) +
        labs(x = "tumor samples",
             y = "mutational signature\nlevel",
             title = glue("Mutational signature levels in {cancer} separated by *KRAS* allele"))
    fn <- as.character(glue("mutsig-memosort_{cancer}_barplot.svg"))
    ggsave_wrapper(p,
                   plot_path(GRAPHS_DIR, fn),
                   width = 10, height = 2)
}

mutsig_noartifact_df %>%
    filter(!all_zeros) %>%
    strip_ras() %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(plot_memosorted_tumor_samples)




#### ---- Simulated artifical signatures ---- ####


d1_function <- function(prob) {
    function(n) { rbinom(n = n, size = 1, prob = prob)}
}

d2_function <- function(mean, sd) {
    function(n) minmax(rnorm(n = n, mean = mean, sd = sd), 0, 1)
}

sampled_signature_groups <- tribble(
    ~ grp, ~ d1_prob, ~ d2_mean, ~ d2_sd,
    "A", 0.1, 0.9, 0.05,
    "B", 0.5, 0.9, 0.05,
    "C", 0.8, 0.5, 0.16,
    "D", 1, 0.5, 0.16
)

sample_signature_levels <- function(d1_prob,
                                    d2_mean, d2_sd,
                                    num_samples,
                                    ...) {
    tibble(ts = seq(1, num_samples)) %>%
        mutate(d1_val = d1_function(d1_prob)(num_samples),
               d2_val = d2_function(d2_mean, d2_sd)(num_samples),
               mut_sig_level = d1_val * d2_val)
}



set.seed(0)
sampled_signature_levels <- sampled_signature_groups
sampled_signature_levels$data <- pmap(
    sampled_signature_groups,
    sample_signature_levels,
    num_samples = 1000
)

sampled_signature_levels %<>% unnest(data)

artificial_mutsig_ggridges <- sampled_signature_levels %>%
    mutate(y_lbl = glue("group: {grp}\nbinomial({d1_prob})\nnormal({d2_mean}, {d2_sd})")) %>%
    ggplot(aes(x = mut_sig_level, y = y_lbl)) +
    geom_density_ridges(aes(color = grp, fill = grp),
                        alpha = 0.5, scale = 1.2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none",
          axis.text.y = element_text(vjust = 0)) +
    labs(x = "simulated mutational signature level",
         y = "distribution density",
         title = "Four different artifical distributions of a mutational signature",
         subtitle = "The distributions are drawn from a mixture of a\nbinomial and normal distribution.")
ggsave_wrapper(artificial_mutsig_ggridges,
               plot_path(GRAPHS_DIR, "artificial-mutsig-ggridges.svg"),
               "small")


set.seed(0)
artificial_mutsig_box <- sampled_signature_levels %>%
    group_by(grp) %>%
    ggplot(aes(x = grp, y = mut_sig_level)) +
    geom_violin(aes(color = grp, fill = grp),
                width = 1, alpha = 0.7,
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(color = "black", fill = "white",
                 width = 0.05,
                 alpha = 0.7, outlier.shape = NA) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none",
          axis.ticks = element_blank()) +
    labs(x = "group",
         y = "artificial mutational signature level")
ggsave_wrapper(artificial_mutsig_box,
               plot_path(GRAPHS_DIR, "artificial-mutsig-boxplots.svg"),
               "small")


sampled_signature_testres <- sampled_signature_levels %>%
    rename(ras_allele = grp,
           contribution = mut_sig_level) %>%
    pairwise_ttest_msigs()


artifical_sig_testres <- sampled_signature_testres %>%
    asymmetrise(group1, group2) %>%
    plot_sig_associations_ggasym(fill_br = diff_med,
                                 fill_br_lbl = "difference\nin medians") +
    coord_fixed()
ggsave_wrapper(artifical_sig_testres,
               plot_path(GRAPHS_DIR, "artifical-sig-testres.svg"),
               "small")
