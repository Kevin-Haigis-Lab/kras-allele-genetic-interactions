# Observed vs. predicted KRAS allele frequency


# A tibble with the number of trinuc. contexts for each KRAS allele mutation.
kras_trinucleotide_contexts_counts <- left_join(
    kras_trinucleotide_contexts, tricontext_counts_df,
    by = "context"
)


# Calculate the frequency of each type of mutation that causes a KRAS allele.
predict_kras_allele_frequency <- function(tib, sequencing_type) {
    if (sequencing_type == "exome") {
        counts_df <- kras_trinucleotide_contexts_counts %>%
            select(-genome_count) %>%
            dplyr::rename(tricontext_count = "exome_count")
    } else if (sequencing_type == "genome") {
        counts_df <- kras_trinucleotide_contexts_counts %>%
            select(-exome_count) %>%
            dplyr::rename(tricontext_count = "genome_count")
    } else {
        stop(glue("'{sequencing_type}' is not a viable option for `sequencing_type`."))
    }

    tib %>%
        count(context, tricontext, name = "mut_count") %>%
        right_join(counts_df, by = c("context", "tricontext")) %>%
        mutate(
            tricontext_mut_count = ifelse(is.na(mut_count), 0, mut_count)
        ) %>%
        select(kras_allele, kras_codon, context, tricontext, tricontext_count,
               tricontext_mut_count)
}


# The number of mutations that mimic that of the KRAS allele.
predicted_kras_allele_frequency <- trinucleotide_mutations_df %>%
    filter(cancer != "SKCM" & !is_hypermutant) %>%
    filter(target %in% c("exome", "genome")) %>%
    dplyr::rename(actual_kras_allele = "kras_allele") %>%
    group_by(cancer, dataset, tumor_sample_barcode, target,
             actual_kras_allele) %>%
    nest() %>%
    ungroup() %>%
    mutate(kras_liklihoods = purrr::map2(data, target, predict_kras_allele_frequency)) %>%
    select(-data) %>%
    unnest(kras_liklihoods) %>%
    group_by(tumor_sample_barcode, dataset, target, cancer, actual_kras_allele,
             kras_allele, kras_codon, context, tricontext_count) %>%
    summarise(
        tricontext = paste(tricontext, collapse = ", "),
        tricontext_mut_count = sum(tricontext_mut_count)
    ) %>%
    ungroup()

ProjectTemplate::cache("predicted_kras_allele_frequency",
                       depends = "trinucleotide_mutations_df")


# For each sample, the liklihood of getting each of the KRAS alleles, given
# that they will get a KRAS mutation at a hotspot.
kras_hostspot_probability <- predicted_kras_allele_frequency %>%
    group_by(tumor_sample_barcode) %>%
    mutate(kras_allele_prob = tricontext_mut_count / tricontext_count,
           kras_allele_prob = kras_allele_prob / sum(kras_allele_prob)) %>%
   ungroup()



#### ---- Plotting ---- ####

# Box plot for distribution of liklihood for each allele in each sample.
predicted_kras_allele_frequency_boxplot <- kras_hostspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    mutate(kras_allele = factor(kras_allele,
                                levels = names(short_allele_pal))) %>%
    ggplot(aes(
        x = kras_allele, y = kras_allele_prob,
        color = kras_allele, fill = kras_allele
    )) +
    facet_wrap(. ~ cancer) +
    geom_boxplot(
        alpha = 0.1,
        outlier.shape = NA
    ) +
    scale_fill_manual(values = short_allele_pal) +
    scale_color_manual(values = short_allele_pal) +
    scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1.0),
        strip.background = element_blank()
    ) +
    labs(
        y = "predicted rate of KRAS allele"
    )

ggsave_wrapper(
    predicted_kras_allele_frequency_boxplot,
    plot_path("50_10_observed-predicted-kras-alleles",
              "predicted_kras_allele_frequency_boxplot.svg"),
    size = "large"
)



# Bar plot for the average predicted frequency of KRAS allele.
# Error bars are SEM.
predict_kras_allele_frequency_barplot1 <- kras_hostspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
        sem_kras_allele_prob = sem(kras_allele_prob),
    ) %>%
    ungroup() %>%
    mutate(
        kras_allele = factor(kras_allele, levels = names(short_allele_pal)),
        upper_line_y = avg_kras_allele_prob + sem_kras_allele_prob,
        lower_line_y = avg_kras_allele_prob - sem_kras_allele_prob,
        lower_line_y = ifelse(lower_line_y < 0, 0, lower_line_y)
    ) %>%
    ggplot(aes(
        x = kras_allele, y = avg_kras_allele_prob,
        color = kras_allele, fill = kras_allele
    )) +
    facet_wrap(. ~ cancer, scales = "free") +
    geom_col(alpha = 0.5) +
    geom_errorbar(
        aes(ymax = upper_line_y, ymin = lower_line_y),
        width = 0.2
    ) +
    scale_fill_manual(values = short_allele_pal) +
    scale_color_manual(values = short_allele_pal) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expand_scale(mult = c(0, 0.02))) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1.0),
        strip.background = element_blank()
    ) +
    labs(
        y = "predicted rate of KRAS allele"
    )
ggsave_wrapper(
    predict_kras_allele_frequency_barplot1,
    plot_path("50_10_observed-predicted-kras-alleles",
              "predict_kras_allele_frequency_barplot1.svg"),
    size = "wide"
)


# The same bar plot as above, but separated by allele to compare across cancer.
predict_kras_allele_frequency_barplot2 <- kras_hostspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
        sem_kras_allele_prob = sem(kras_allele_prob),
    ) %>%
    ungroup() %>%
    mutate(
        kras_allele = factor(kras_allele, levels = names(short_allele_pal)),
        upper_line_y = avg_kras_allele_prob + sem_kras_allele_prob,
        lower_line_y = avg_kras_allele_prob - sem_kras_allele_prob,
        lower_line_y = ifelse(lower_line_y < 0, 0, lower_line_y)
    ) %>%
    ggplot(aes(
        x = cancer, y = avg_kras_allele_prob,
        color = cancer, fill = cancer
    )) +
    facet_wrap(~ kras_allele, scales = "free") +
    geom_col(alpha = 0.5) +
    geom_errorbar(
        aes(ymax = upper_line_y, ymin = lower_line_y),
        width = 0.2
    ) +
    scale_fill_manual(values = cancer_palette) +
    scale_color_manual(values = cancer_palette) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expand_scale(mult = c(0, 0.02))) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1.0),
        strip.background = element_blank()
    ) +
    labs(
        y = "predicted rate of KRAS allele"
    )

ggsave_wrapper(
    predict_kras_allele_frequency_barplot2,
    plot_path("50_10_observed-predicted-kras-alleles",
              "predict_kras_allele_frequency_barplot2.svg"),
    size = "large"
)


# Make a single scatter plot for an individual cancer.
# This was made to be called from `obs_v_pred_scatter_plot` for each cancer.
individual_obs_v_pred_scatter_plot <- function(plot_df, cancer, stats = FALSE) {

    axis_lim <- max(c(plot_df$real_kras_allele_frequency,
                      plot_df$avg_kras_allele_prob))

    p <- ggplot(plot_df,
                aes(x = real_kras_allele_frequency,
                    y = avg_kras_allele_prob))

    if (stats) {
        p <- p + geom_point(aes(color = kras_allele,
                                size = log_p_value,
                                shape = point_shape),
                            alpha = 0.8) +
            scale_shape_identity() +
            scale_size_continuous(
                range = c(1, 3),
                guide = guide_legend(
                    title.position = "left",
                    keywidth = unit(4, "mm"),
                    label.position = "top",
                    nrow = 1,
                    order = 1
            ))
    } else {
        p <- p + geom_point(aes(color = kras_allele),
                            size = 1, alpha = 0.7)
    }

    p <- p +
        geom_abline(intercept = 0, slope = 1) +
        ggrepel::geom_text_repel(
            aes(label = kras_allele),
            color = "grey25",
            family = "Arial",
            size = 2.5,
            force = 0.1,
            segment.alpha = 0.5,
            segment.size = 0.1,
            seed = 0
        ) +
        scale_color_manual(
            values = short_allele_pal,
            guide = FALSE
        ) +
        scale_x_continuous(
            limits = c(0, axis_lim),
            expand = expand_scale(mult = c(0, 0.02))
        ) +
        scale_y_continuous(
            limits = c(0, axis_lim),
            expand = expand_scale(mult = c(0, 0.02))
        ) +
        coord_fixed() +
        theme_bw(base_size = 8, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.title.align = 0.5,
            legend.spacing = unit(2, "mm"),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2, "mm"),
            strip.background = element_blank()
        ) +
        labs(
            title = cancer,
            color = "KRAS allele",
            size = "-log10( p-value )",
            x = "observed KRAS allele frequency",
            y = "predicted KRAS allele frequency"
        )
    return(p)
}


# Scatter plot of observed vs. predicted KRAS allele frequency.
obs_v_pred_scatter_plot <- function(real_tib, pred_tib,
                                    save_template,
                                    stats_tib = NULL, p_value_cut = 0.05,
                                    save_size = "small") {
    pdata <- left_join(real_tib, pred_tib,
                        by = c("cancer", "kras_allele")) %>%
        filter(kras_allele %in% names(short_allele_pal)) %>%
        filter(!is.na(avg_kras_allele_prob)) %>%
        filter(cancer != "SKCM")

    if (!is.null(stats_tib)) {
        pdata <- stats_tib %>%
            select(cancer, kras_allele, p_value) %>%
            right_join(pdata, by = c("cancer", "kras_allele")) %>%
            mutate(log_p_value = -log10(p_value),
                   point_shape = ifelse(p_value < 0.05, 16, 17))
    }

    for (CANCER in pdata$cancer) {
        p <- pdata %>%
            filter(cancer == !!CANCER) %>%
            individual_obs_v_pred_scatter_plot(cancer = CANCER,
                                               stats = !is.null(stats_tib))
        ggsave_wrapper(
            p,
            plot_path("50_10_observed-predicted-kras-alleles",
                      glue(save_template)),
            save_size
        )
    }
}




# Scatter plot of all alleles observed vs. predicted KRAS allele frequency.
real_af <- cancer_muts_df %>%
    filter(ras == "KRAS") %>%
    select(cancer, tumor_sample_barcode, ras_allele) %>%
    group_by(cancer, tumor_sample_barcode, ras_allele) %>%
    slice(1) %>%
    ungroup() %>%
    unique() %>%
    count(cancer, ras_allele, name = "ras_allele_count") %>%
    group_by(cancer) %>%
    mutate(real_kras_allele_frequency = ras_allele_count / sum(ras_allele_count)) %>%
    ungroup() %>%
    dplyr::rename(kras_allele = "ras_allele") %>%
    mutate(kras_allele = str_remove(kras_allele, "KRAS_"))

predicted_af <- kras_hostspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
        sem_kras_allele_prob = sem(kras_allele_prob),
    ) %>%
    ungroup()

obs_v_pred_scatter_plot(real_af, predicted_af,
                        save_template = "obs_pred_plot_{CANCER}.svg")



# Scatter plot of G12 alleles observed vs. predicted KRAS allele frequency.
predicted_af_g12 <- predicted_kras_allele_frequency %>%
    filter(kras_codon == 12) %>%
    group_by(tumor_sample_barcode) %>%
    mutate(
        kras_allele_prob = tricontext_mut_count / tricontext_count,
        kras_allele_prob = kras_allele_prob / sum(kras_allele_prob)
    ) %>%
    filter(!is.na(kras_allele_prob)) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
        sem_kras_allele_prob = sem(kras_allele_prob),
    ) %>%
    ungroup()

obs_v_pred_scatter_plot(real_af, predicted_af_g12,
                        save_template = "obs_pred_plot_g12_{CANCER}.svg")




#### ---- Statisics ---- ####

cancer_sample_count_df <- cancer_muts_df %>%
    filter(ras_allele != "WT") %>%
    group_by(cancer, tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    unique() %>%
    count(cancer, name = "num_cancer_samples")

# Binomial test for the actual vs. predicted rates of KRAS alleles.
kras_allele_binomial_test <- function(df) {
    b <- binom.test(x = df$ras_allele_count[[1]],
                    n = df$num_cancer_samples[[1]],
                    p = df$avg_kras_allele_prob[[1]])
    return(tidy(b))
}


kras_allele_freq_stats <- left_join(
        real_af, predicted_af, by = c("cancer", "kras_allele")
    ) %>%
    filter(kras_allele %in% names(short_allele_pal)) %>%
    filter(cancer != "SKCM") %>%
    left_join(cancer_sample_count_df, by = "cancer") %>%
    group_by(cancer, kras_allele) %>%
    nest() %>%
    ungroup() %>%
    mutate(binom_res = purrr::map(data, kras_allele_binomial_test)) %>%
    unnest(c(data, binom_res)) %>%
    janitor::clean_names()

write_tsv(
    kras_allele_freq_stats,
    file.path("tables",
              "50_10_observed-predicted-kras-alleles",
              "kras_allele_freq_stats.tsv")
)

kras_allele_freq_stats %>%
    filter(p_value >= 0.05) %>%
    select(cancer, kras_allele,
           real_kras_allele_frequency, avg_kras_allele_prob,
           p_value)

obs_v_pred_scatter_plot(real_af, predicted_af,
                        stats_tib = kras_allele_freq_stats,
                        save_template = "obs_pred_plot_stats_{CANCER}.svg")




# Just G12 mutations.
kras_g12_freq_stats <- left_join(
        real_af, predicted_af_g12, by = c("cancer", "kras_allele")
    ) %>%
    filter(kras_allele %in% names(short_allele_pal) & !is.na(avg_kras_allele_prob)) %>%
    filter(cancer != "SKCM") %>%
    left_join(cancer_sample_count_df, by = "cancer") %>%
    group_by(cancer, kras_allele) %>%
    nest() %>%
    ungroup() %>%
    mutate(binom_res = purrr::map(data, kras_allele_binomial_test)) %>%
    unnest(c(data, binom_res)) %>%
    janitor::clean_names()

write_tsv(
    kras_g12_freq_stats,
    file.path("tables",
              "50_10_observed-predicted-kras-alleles",
              "kras_g12_freq_stats.tsv")
)

kras_g12_freq_stats %>%
    filter(p_value >= 0.05) %>%
    select(cancer, kras_allele,
           real_kras_allele_frequency, avg_kras_allele_prob,
           p_value)

obs_v_pred_scatter_plot(real_af, predicted_af_g12,
                        stats_tib = kras_g12_freq_stats,
                        save_template = "obs_pred_plot_g12_stats_{CANCER}.svg")
