
GRAPHS_DIR <- "90_03_mutation-burden-distribution"
reset_graph_directory(GRAPHS_DIR)

theme_set(theme_bw(base_size = 7, base_family = "Arial"))


#### ---- Mut. burden distribution ---- ####
# Plots showing the distribution of the mutations per sample.

# Plot the distribution of mutational burden per sample.
# (Also saves to SVG and proto for figures.)

plot_distribution_of_mutation_count <- function(cancer,
                                                data,
                                                save_name,
                                                hypermutant_rug = FALSE) {
    p <- data %>%
        mutate(tumor_sample_barcode = fct_reorder(tumor_sample_barcode, n),
               dataset_label = paste0(dataset, "\n", target)) %>%
        ggplot(
            aes(x = tumor_sample_barcode,
                y = log10(n))
        ) +
        facet_grid(~ dataset_label, scales = "free_x") +
        geom_point(
            aes(color = log10(n)),
            size = 0.1
        ) +
        scale_color_viridis_c(
            guide = FALSE
        ) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.title.y = element_markdown(),
            panel.grid.major.x = element_blank(),
            strip.background = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        labs(
            title = as.character(
                glue("Distribution of mutations in {cancer} samples")
            ),
            x = "tumor samples",
            y = "num. mutations (*log*<sub>10</sub>-transformed)"
        )

    if (hypermutant_rug) {
        p <- p +
            geom_rug(
                aes(alpha = is_hypermutant),
                sides = "b",
                color = "black"
            ) +
            scale_alpha_manual(
                values = c("TRUE" = 0.2, "FALSE" = 0.0),
                guide = FALSE
            )
    }
    save_path <- plot_path(GRAPHS_DIR, save_name)
    ggsave_wrapper(p, save_path, width = 8, height = 2.5)

    # Save for use in a figure.
    saveFigRds(p, save_path)
}

# Distribution of all mutations per sample.
cancer_muts_df %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    group_by(cancer, tumor_sample_barcode, dataset, target,
             is_hypermutant, ras_allele) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_all_muts_distribution.svg")) %>%
    pwalk(plot_distribution_of_mutation_count)

# Distribution of coding mutations per sample.
cancer_coding_muts_df %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    group_by(cancer, tumor_sample_barcode, dataset, target,
             is_hypermutant, ras_allele) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(
        save_name = paste0(cancer, "_coding_muts_distribution.svg"),
        hypermutant_rug = cancer == "COAD"
    ) %>%
    pwalk(plot_distribution_of_mutation_count)



#### ---- Distribution of mutation types per data source ---- ####

plot_distribution_of_mutation_type <- function(cancer, data, save_name) {
    # Get the ordering for the mutation types.
    mutation_type_hr_order <- data %>%
        group_by(mutation_type_hr) %>%
        summarise(tot = sum(n)) %>%
        ungroup() %>%
        arrange(tot) %>%
        pull(mutation_type_hr) %>%
        unlist()

    p <- data %>%
        mutate(
            dataset_label = paste0(dataset, "\n", target),
            mutation_type_hr = factor(mutation_type_hr,
                                      levels = mutation_type_hr_order)
        ) %>%
        ggplot(aes(x = dataset_label, y = n)) +
        geom_col(
            aes(fill = mutation_type_hr),
            position = "fill"
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(
            values = mutation_pal,
            guide = guide_legend(
                nrow = 2
            )
        ) +
        theme_bw(base_size = 8, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 40, hjust = 1.0),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        labs(
            y = "fraction of mutations",
            fill = "mutation type"
        )

    save_path <- plot_path("90_03_mutation-burden-distribution", save_name)
    ggsave_wrapper(p, save_path, width = 6, height = 4)

    # Save the ggplot objects for figures.
    saveFigRds(p, save_path)

    return(p)
}


# Distribution of coding mutations per sample.
mut_type_bars <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM") %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    group_by(cancer, dataset, target, mutation_type_hr) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_mutation_types.svg")) %>%
    pwalk(plot_distribution_of_mutation_type)



#### ---- VAF distribution ---- ####

# A function to make the labels of a ggplot a bit prettier.
pretty_rounding <- function(x) {
    case_when(x == 0 ~ "0",
              x == 1 ~ "1",
              TRUE ~ as.character(round(x, 2)))
}


# Density plot of the distribution of VAFs.
vaf_distribution_density <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    ggplot(aes(x = VAF)) +
    facet_wrap(~ cancer, scales = "free_y") +
    geom_density(aes(color = cancer)) +
    scale_color_manual(values = cancer_palette, guide = FALSE) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(
        limits = c(0, 1),
        expand = c(0, 0),
        labels = pretty_rounding
    ) +
    theme_bw(base_size = 7, base_family = "arial") +
    theme(
        strip.background = element_blank()
    )
ggsave_wrapper(
    vaf_distribution_density,
    plot_path(GRAPHS_DIR, "vaf_distribution_density.svg"),
    "small"
)



plot_mutfreq_by_avgvaf_scatter <- function(df) {
    df %>%
        ggplot(aes(x = avg_vaf, y = mut_freq)) +
        facet_wrap(~ cancer, scales = "free") +
        geom_point(size = 0.3, alpha = 0.7, color = "grey25") +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 1.2, family = "arial", color = "blue",
            segment.color = "grey40", segment.size = 0.2, segment.alpha = 0.5,
            seed = 0,
        ) +
        scale_x_continuous(
            limits = c(0, 1),
            expand = c(0, 0),
            labels = pretty_rounding
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.02)),
            labels = pretty_rounding
        ) +
        theme_bw(base_size = 7, base_family = "arial") +
        theme(
            strip.background = element_blank()
        ) +
        labs(
            x = "average VAF",
            y = "mutation frequency"
        )
}

# Average VAF of mutations in a gene vs. mutational frequency of the gene.
vaf_mutfreq_scatter <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    group_by(cancer) %>%
    mutate(num_cancer_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer, hugo_symbol) %>%
    summarise(
        avg_vaf = mean(VAF),
        mut_freq = n_distinct(tumor_sample_barcode) / unique(num_cancer_samples)
    ) %>%
    group_by(cancer) %>%
    arrange(-mut_freq, -avg_vaf) %>%
    mutate(label = ifelse(1:n() < 20, hugo_symbol, NA)) %>%
    ungroup() %>%
    plot_mutfreq_by_avgvaf_scatter()

ggsave_wrapper(
    vaf_mutfreq_scatter,
    plot_path(GRAPHS_DIR, "vaf_mutfreq_scatter.jpeg"),
    "medium"
)


# A simple data frame of the genes with genetic interactions per cancer type.
genetic_interaction_genes <- genetic_interaction_df %>%
    select(cancer, hugo_symbol) %>%
    unique() %>%
    bind_rows(
        tibble(cancer = names(cancer_palette), hugo_symbol = "KRAS")
    )


# Average VAF of mutations in a gene vs. mutational frequency of the gene, but
# only for genes with a comutation interaction
vaf_mutfreq_comuts_scatter <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    inner_join(genetic_interaction_genes, by = c("cancer", "hugo_symbol")) %>%
    group_by(cancer) %>%
    mutate(num_cancer_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer, hugo_symbol) %>%
    summarise(
        avg_vaf = mean(VAF),
        mut_freq = n_distinct(tumor_sample_barcode) / unique(num_cancer_samples)
    ) %>%
    group_by(cancer) %>%
    arrange(-mut_freq, -avg_vaf) %>%
    mutate(label = ifelse(1:n() < 20, hugo_symbol, NA)) %>%
    ungroup() %>%
    plot_mutfreq_by_avgvaf_scatter()

ggsave_wrapper(
    vaf_mutfreq_comuts_scatter,
    plot_path(GRAPHS_DIR, "vaf_mutfreq_comuts_scatter.svg"),
    "medium"
)


# Summary statistics on VAF values per cancer type.
cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    group_by(cancer) %>%
    summarise(
        min_vaf = min(VAF),
        q25_vaf = quantile(VAF, 0.25),
        avg_vaf = mean(VAF),
        mid_vaf = median(VAF),
        q75_vaf = quantile(VAF, 0.75),
        max_vaf = max(VAF),
        stddev_vaf = sd(VAF)
    ) %>%
    ungroup() %>%
    knitr::kable(digits = 3)


#### ---- VAF figures for resubmission ---- ####

distribution_kras_vaf <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    inner_join(tcga_purity_ploidy, by = "tumor_sample_barcode") %>%
    filter(!is.na(purity)) %>%
    filter(ras_allele != "WT") %>%
    filter(hugo_symbol == "KRAS") %>%
    mutate(adjusted_vaf = VAF / cancer_dna_fraction,
           adjusted_vaf = map_dbl(adjusted_vaf, ~ min(.x, 1))) %>%
    filter(hugo_symbol == "KRAS") %>%
    ggplot(aes(adjusted_vaf)) +
    facet_wrap(~ cancer, nrow = 1, scales = "free_y") +
    geom_histogram(aes(y = ..density.., fill = cancer, color = cancer),
                   breaks = seq(0, 1, 0.05),
                   size = 0.2, alpha = 0.5) +
    geom_density(aes(color = cancer), size = 1, fill = NA) +
    scale_color_manual(values = cancer_palette, drop = TRUE) +
    scale_fill_manual(values = cancer_palette, drop = TRUE) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          panel.spacing.x = unit(3, "mm"),
          axis.title.x = element_markdown()) +
    labs(x = "*KRAS* mutation VAF (adjusted for tumor sample purity)",
         y = "density")
ggsave_wrapper(distribution_kras_vaf,
               plot_path(GRAPHS_DIR, "kras-adj-vaf-distribution.svg"),
               "wide")
saveFigRds(distribution_kras_vaf, "kras-adj-vaf-distribution.svg")



comut_genes_vaf_dist <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    inner_join(tcga_purity_ploidy, by = "tumor_sample_barcode") %>%
    filter(!is.na(purity)) %>%
    filter(hugo_symbol != "KRAS") %>%
    mutate(adjusted_vaf = VAF / cancer_dna_fraction,
           adjusted_vaf = map_dbl(adjusted_vaf, ~ min(.x, 1))) %>%
    left_join(
        genetic_interaction_df %>%
            select(hugo_symbol, cancer) %>%
            add_column(is_comut = TRUE) %>%
            distinct(),
        by = c("cancer", "hugo_symbol")
    ) %>%
    mutate(
        is_comut = ifelse(is.na(is_comut), FALSE, is_comut),
        is_comut = ifelse(is_comut,
                          "of comutating gene",
                          "of non-comutating gene")
    ) %>%
    ggplot(aes(adjusted_vaf)) +
    facet_wrap(~ cancer, nrow = 1, scales = "free_y") +
    geom_density(aes(color = is_comut, fill = is_comut),
                 size = 0.9, alpha = 0.2) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          panel.spacing.x = unit(3, "mm"),
          axis.title.x = element_markdown()) +
    labs(x = "mutation VAF (adjusted for tumor sample purity)",
         y = "density")
ggsave_wrapper(comut_genes_vaf_dist,
               plot_path(GRAPHS_DIR, "comutation-genes-vaf-dist.svg"),
               "wide")
saveFigRds(comut_genes_vaf_dist, "comutation-genes-vaf-dist")


cancer_muts_adjvaf <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    inner_join(tcga_purity_ploidy, by = "tumor_sample_barcode") %>%
    filter(!is.na(purity)) %>%
    mutate(adjusted_vaf = VAF / cancer_dna_fraction,
           adjusted_vaf = map_dbl(adjusted_vaf, ~ min(.x, 1)))


kras_mut_adjvaf <- cancer_muts_adjvaf %>%
    filter(ras_allele != "WT") %>%
    filter(hugo_symbol == "KRAS") %>%
    select(tumor_sample_barcode, cancer, kras_adj_vaf = adjusted_vaf)


diff_adjusted_vaf_dist <- cancer_muts_adjvaf %>%
    filter(ras_allele != "WT") %>%
    filter(hugo_symbol != "KRAS") %>%
    inner_join(kras_mut_adjvaf, by = c("tumor_sample_barcode", "cancer")) %>%
    inner_join(
        genetic_interaction_df %>%
            select(hugo_symbol, cancer) %>%
            distinct(),
        by = c("cancer", "hugo_symbol")
    ) %>%
    mutate(diff_adj_vaf = adjusted_vaf - kras_adj_vaf) %>%
    ggplot(aes(diff_adj_vaf)) +
    facet_wrap(~ cancer, nrow = 1, scales = "free_y") +
    geom_density(aes(color = cancer, fill = cancer),
                 size = 0.9, alpha = 0.2) +
    scale_color_manual(values = cancer_palette, drop = TRUE) +
    scale_fill_manual(values = cancer_palette, drop = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          panel.spacing.x = unit(3, "mm"),
          axis.title.x = element_markdown()) +
    labs(x = "difference in VAF of *KRAS* and comutated gene mutations",
         y = "density")
ggsave_wrapper(diff_adjusted_vaf_dist,
               plot_path(GRAPHS_DIR, "kras-comuts-vaf-difference.svg"),
               "wide")



# cancer_muts_adjvaf %>%
#     filter(!is_hypermutant) %>%
#     filter(cancer == "PAAD") %>%
#     saveRDS("PAAD_cancer_muts_adjvaf.rds")


# cancer_full_muts_df %>%
#     filter(!is_hypermutant) %>%
#     filter(str_detect(mutation_type, "synon|silent")) %>%
#     filter(!is.na(VAF)) %>%
#     filter(cancer == "COAD") %>%
#     saveRDS("cancer_silent_muts_df.rds")


filter_cancers <- function(df) {
    df %>%
        filter(cancer %in% c("COAD", "LUAD", "PAAD"))
}


density_margins_plot <- function(df, x, y, x_label, y_label) {
    margin_plot <- function(df, x) {
        df %>%
            ggplot(aes({{ x }})) +
            geom_density(color = "dodgerblue", fill = "dodgerblue",
                         size = 1.2, alpha = 0.4) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
    }

    p1 <- df %>%
        ggplot(aes({{ x }}, {{ y }})) +
        geom_jitter(alpha = 0.1, size = 0.05, width = 0.01, height = 0.01) +
        # geom_density2d(color = "tomato", alpha = 0.9, size = 1.2,
        #                binwidth = 0.05) +
        labs(x = x_label,
             y = y_label)

    p2 <- margin_plot(df, {{ x }})+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
    p3 <- margin_plot(df, {{ y }}) +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank()) +
        coord_flip()


    p_design <- c(
        area(t = 1, l = 1, b = 2, r = 5),
        area(t = 3, l = 1, b = 7, r = 5),
        area(t = 3, l = 6, b = 7, r = 7)
    )
    patch <- p2 + p1 + p3 +
        plot_layout(design = p_design) &
        theme(axis.ticks = element_blank())
    return(patch)
}

density_margins_per_cancer <- function(cancer, data, x, y, x_label, y_label,
                                       plt_name_glue) {
    p <- density_margins_plot(
        data,
        x = {{ x }},
        y = {{ y }},
        x_label = x_label,
        y_label = y_label
    ) +
        plot_annotation(title = cancer,
                        theme = theme(plot.title = element_text(hjust = 0.5)))

    plt_name <- as.character(glue(plt_name_glue))
    ggsave_wrapper(p,
                   plot_path(GRAPHS_DIR, plt_name),
                   "medium")
}

cancer_muts_adjvaf %>%
    filter_cancers() %>%
    group_by(cancer) %>%
    sample_n(2e4) %>%
    nest() %>%
    pwalk(
        density_margins_per_cancer,
        x = purity,
        y = VAF,
        x_label = "tumor sample purity",
        y_label = "VAF (unadjusted)",
        plt_name_glue = "vaf-density2d_{cancer}.svg"
    )

cancer_muts_adjvaf %>%
    filter_cancers() %>%
    group_by(cancer) %>%
    sample_n(2e4) %>%
    nest() %>%
    pwalk(
        density_margins_per_cancer,
        x = purity,
        y = adjusted_vaf,
        x_label = "tumor sample purity",
        y_label = "VAF (adjusted for tumor purity)",
        plt_name_glue = "vaf-adjusted-density2d_{cancer}.svg"
    )




subtitle <- "The dashed line represents where the adjustment for tumor purity had no effect.
Each point represents the mean value for a single tumor sample."

vaf_v_adjvaf_scatter <- cancer_muts_adjvaf %>%
    group_by(cancer, tumor_sample_barcode) %>%
    summarise(avg_VAF = mean(VAF),
              avg_adjusted_vaf = mean(adjusted_vaf),
              purity = mean(purity)) %>%
    ggplot(aes(avg_VAF, avg_adjusted_vaf)) +
    geom_abline(slope = 1, intercept = 0, lty = 2, color = "grey50") +
    geom_jitter(aes(size = purity, color = cancer),
                alpha = 0.5, height = 0.02, width = 0.02) +
    scale_size_continuous(range = c(1, 5)) +
    scale_color_manual(values = cancer_palette, drop = TRUE) +
    labs(x = "VAF",
         y = "VAF adjusted for tumor sample purity",
         size = "tumor sample purity",
         color = "cancer type",
         title = "Effects of adjusting VAF for tumor sample purity",
         subtitle = subtitle)


ggsave_wrapper(vaf_v_adjvaf_scatter,
               plot_path(GRAPHS_DIR, "vaf-v-adjvaf-scatter.svg"),
               "medium")


jitter_pos <- position_jitter(width = 0.02, height = 0, seed = 0)
subtitle <- "Each point represents the median (± the std. dev.) VAF of an individual tumor sample."
purity_vs_vaf_scatter <- cancer_muts_adjvaf %>%
    group_by(cancer, tumor_sample_barcode) %>%
    summarise(avg_VAF = median(VAF),
              sd_VAF = sd(VAF),
              purity = mean(purity)) %>%
    ungroup() %>%
    mutate(upper_VAF = minmax(avg_VAF + sd_VAF, 0, 1),
           lower_VAF = minmax(avg_VAF - sd_VAF, 0, 1)) %>%
    ggplot(aes(x = purity, y = avg_VAF)) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_linerange(aes(ymin = lower_VAF, ymax = upper_VAF),
                   alpha = 0.3, position = jitter_pos) +
    geom_point(position = jitter_pos) +
    geom_smooth(method = "lm", formula = "y ~ x") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    labs(x = "tumor sample purity",
         y = "mutation VAF (unadjusted for tumor sample purity)",
         title = "Correlation between tumor sample purity and VAF values",
         subtitle = subtitle)
ggsave_wrapper(purity_vs_vaf_scatter,
               plot_path(GRAPHS_DIR, "purity-vs-vaf-scatter.svg"),
               "large")
