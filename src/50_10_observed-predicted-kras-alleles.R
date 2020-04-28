# Observed vs. predicted KRAS allele frequency

GRAPHS_DIR <- "50_10_observed-predicted-kras-alleles"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


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
        stop(glue(
            "'{sequencing_type}' is not a viable option for `sequencing_type`."
        ))
    }

    tib %>%
        count(context, tricontext, name = "mut_count") %>%
        right_join(counts_df, by = c("context", "tricontext")) %>%
        mutate(
            total_num_mutations = sum(mut_count, na.rm = TRUE),
            tricontext_mut_count = ifelse(is.na(mut_count), 0, mut_count)
        ) %>%
        select(kras_allele, kras_codon, context, tricontext, tricontext_count,
               tricontext_mut_count, total_num_mutations)
}


# The number of mutations that mimic that of the KRAS allele.
ProjectTemplate::cache("predicted_kras_allele_frequency",
                       depends = "trinucleotide_mutations_df",
{
    predicted_kras_allele_frequency <- trinucleotide_mutations_df %>%
        filter(cancer != "SKCM" & !is_hypermutant) %>%
        filter(target %in% c("exome", "genome")) %>%
        filter(hugo_symbol != "KRAS") %>%
        dplyr::rename(actual_kras_allele = "kras_allele") %>%
        group_by(cancer, dataset, tumor_sample_barcode, target,
                 actual_kras_allele) %>%
        nest() %>%
        ungroup() %>%
        mutate(kras_liklihoods = purrr::map2(
            data, target, predict_kras_allele_frequency
        )) %>%
        select(-data) %>%
        unnest(kras_liklihoods) %>%
        group_by(tumor_sample_barcode, dataset, target, cancer,
                 actual_kras_allele, kras_allele, kras_codon,
                 context, tricontext_count, total_num_mutations) %>%
        summarise(
            tricontext = paste(tricontext, collapse = ", "),
            tricontext_mut_count = sum(tricontext_mut_count)
        ) %>%
        ungroup()
    return(predicted_kras_allele_frequency)
})


alleles_frequency_per_cancer_df <- predicted_kras_allele_frequency %>%
    select(tumor_sample_barcode, cancer, actual_kras_allele) %>%
    unique() %>%
    count(cancer, actual_kras_allele) %>%
    dplyr::rename(kras_allele = actual_kras_allele)


# For each sample, the liklihood of getting each of the KRAS alleles, given
# that they will get a KRAS mutation at a hotspot.
remove_alleles_with_low_frequency <- function(tib, min_num = 15) {

    if (!"kras_allele" %in% colnames(tib)) {
        stop("There is no column 'kras_allele' in `tib`.")
    }

    left_join(
        tib, alleles_frequency_per_cancer_df,
        by = c("cancer", "kras_allele")
    ) %>%
        filter(n >= !!min_num) %>%
        select(-n)
}

kras_hotspot_probability <- predicted_kras_allele_frequency %>%
    remove_alleles_with_low_frequency() %>%
    group_by(tumor_sample_barcode) %>%
    mutate(
        kras_allele_prob = tricontext_mut_count / tricontext_count,
        kras_allele_prob = kras_allele_prob / sum(kras_allele_prob),
    ) %>%
    ungroup()



#### ---- Bootstrap 95% CI for predictions ---- ####

# A function to be passed as `statistic` to `boot::boot()`.
# `df` is the original data frame with one row per data point.
# `index` is the index to use for the sampling process.
#
# This function uses the already calculated frequency of each KRAS allele-type
#   mutation in each sample and just re-calculates the average prediction across
#   all samples.
boot_predict_kras_allele_frequency <- function(df, index) {
    mod_df <- df[index, ]
    mod_df$tumor_sample_barcode <- paste0(mod_df$tumor_sample_barcode, "_", index)
    res <- mod_df %>%
        unnest(data) %>%
        group_by(tumor_sample_barcode) %>%
        mutate(
            kras_allele_prob = tricontext_mut_count / tricontext_count,
            kras_allele_prob = kras_allele_prob / sum(kras_allele_prob, na.rm = TRUE),
        ) %>%
        ungroup() %>%
        group_by(kras_allele) %>%
        summarise(avg_kras_allele_prob = mean(kras_allele_prob,
                                              na.rm = TRUE)) %>%
        ungroup() %>%
        select(kras_allele, avg_kras_allele_prob) %>%
        deframe()
    return(res)
}

# Run the bootstrapping process for a cancer.
# `data` should be a data frame with one row per tumor sample.
bootstrap_allele_confidence_intervals <- function(cancer, data, R = 1e3) {
    grouped_data <- data %>%
        group_by(tumor_sample_barcode) %>%
        nest() %>%
        ungroup()

    allele_bs <- boot::boot(
        data = grouped_data,
        statistic = boot_predict_kras_allele_frequency,
        R = R
    )

    result_tib <- tibble(
            kras_allele = names(allele_bs$t0)
        ) %>%
        mutate(
            index = seq(1, n(), 1),
            ci_obj = purrr::map(index, ~ boot::boot.ci(allele_bs,
                                                       type = "perc",
                                                       index = .x)),
            ci_upper = purrr::map_dbl(ci_obj, ~ .x$percent[[4]]),
            ci_lower = purrr::map_dbl(ci_obj, ~ .x$percent[[5]])
        )
    return(result_tib)
}

# Run the bootstrapping process to estimate the 95% CI for the KRAS allele
# frequency prediction.
ProjectTemplate::cache("kras_allele_freq_bootstrap_ci",
                       depends = "predicted_kras_allele_frequency",
{
    kras_allele_freq_bootstrap_ci <- predicted_kras_allele_frequency %>%
        remove_alleles_with_low_frequency() %>%
        group_by(cancer) %>%
        nest() %>%
        ungroup() %>%
        mutate(
            bs_results = purrr::map2(cancer, data,
                                     bootstrap_allele_confidence_intervals)
        ) %>%
        select(-data) %>%
        unnest(bs_results)
    return(kras_allele_freq_bootstrap_ci)
})




#### ---- Plotting ---- ####

# Box plot for distribution of likelihood for each allele in each sample.
predicted_kras_allele_frequency_boxplot <- kras_hotspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    mutate(kras_allele = factor_alleles(kras_allele)) %>%
    ggplot(aes(
        x = kras_allele, y = kras_allele_prob,
        color = kras_allele, fill = kras_allele
    )) +
    facet_wrap(. ~ cancer, scales = "free_x") +
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
    plot_path(GRAPHS_DIR,
              "predicted_kras_allele_frequency_boxplot.svg"),
    size = "large"
)



# Bar plot for the average predicted frequency of KRAS allele.
# Error bars are SEM.
barplot_df <- kras_hotspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    mutate(
        kras_allele_prob_weighted = kras_allele_prob * log2(total_num_mutations)
    ) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
        avg_kras_allele_prob_weighted = mean(kras_allele_prob_weighted)
    ) %>%
    group_by(cancer) %>%
    mutate(
        avg_kras_allele_prob = avg_kras_allele_prob / sum(avg_kras_allele_prob),
        avg_kras_allele_prob_weighted = avg_kras_allele_prob_weighted / sum(avg_kras_allele_prob_weighted)
    ) %>%
    ungroup() %>%
    left_join(kras_allele_freq_bootstrap_ci,
              by = c("cancer", "kras_allele")) %>%
    mutate(
        ci_upper_avg_kras_allele_prob = avg_kras_allele_prob + ci_upper,
        ci_lower_avg_kras_allele_prob = avg_kras_allele_prob - ci_lower
    )


predict_kras_allele_frequency_barplot1 <- barplot_df %>%
    mutate(
        kras_allele = factor_alleles(kras_allele)
    ) %>%
    ggplot(aes(
        x = kras_allele, y = avg_kras_allele_prob,
        color = kras_allele, fill = kras_allele
    )) +
    facet_wrap(. ~ cancer, scales = "free") +
    geom_col(alpha = 0.5) +
    geom_errorbar(
        aes(ymin = ci_lower_avg_kras_allele_prob,
            ymax = ci_upper_avg_kras_allele_prob),
        width = 0.2
    ) +
    scale_fill_manual(values = short_allele_pal) +
    scale_color_manual(values = short_allele_pal) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.02))) +
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
    plot_path(GRAPHS_DIR, "predict_kras_allele_frequency_barplot1.svg"),
    size = "wide"
)

predict_kras_allele_frequency_weighted_barplot1 <- barplot_df %>%
    mutate(
        kras_allele = factor_alleles(kras_allele)
    ) %>%
    ggplot(aes(
        x = kras_allele, y = avg_kras_allele_prob_weighted,
        color = kras_allele, fill = kras_allele
    )) +
    facet_wrap(. ~ cancer, scales = "free") +
    geom_col(alpha = 0.5) +
    scale_fill_manual(values = short_allele_pal) +
    scale_color_manual(values = short_allele_pal) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.02))) +
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
    predict_kras_allele_frequency_weighted_barplot1,
    plot_path(GRAPHS_DIR,
              "predict_kras_allele_frequency_weighted_barplot1.svg"),
    size = "wide"
)


single_kras_allele_freq_barplot <- function(cancer, data,
                                            y_max = NA,
                                            with_errorbars = TRUE) {
    p <- data %>%
        mutate(
            kras_allele = fct_reorder(kras_allele, avg_kras_allele_prob)
        ) %>%
        ggplot(aes(
            x = kras_allele, y = avg_kras_allele_prob,
            color = kras_allele, fill = kras_allele
        )) +
        geom_col(alpha = 0.5)

    if (with_errorbars) {
        p <- p + geom_errorbar(
                               aes(ymin = ci_lower_avg_kras_allele_prob,
                                   ymax = ci_upper_avg_kras_allele_prob),
                               width = 0.2
                               )
    }

    p <- p +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        scale_y_continuous(limits = c(0, y_max),
                           expand = expansion(mult = c(0, 0.02))) +
        theme_bw(base_size = 8, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1.0),
            strip.background = element_blank()
        ) +
        labs(title = cancer)
    return(p)
}

# Bar plots of the un-weighted predicted frequencies.
predict_kras_allele_frequency_barplot2_df <- barplot_df %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(barplot = purrr::map2(
        cancer, data, single_kras_allele_freq_barplot,
        y_max = max(barplot_df$ci_upper_avg_kras_allele_prob)
    ))

predict_kras_allele_frequency_barplot2 <- patchwork::wrap_plots(
    predict_kras_allele_frequency_barplot2_df$barplot,
    nrow = 2
)

ggsave_wrapper(
    predict_kras_allele_frequency_barplot2,
    plot_path(GRAPHS_DIR, "predict_kras_allele_frequency_barplot2.svg"),
    "wide"
)


# Use the weighted prediction.
predict_kras_allele_frequency_weighted_barplot2_df <- barplot_df %>%
    mutate(avg_kras_allele_prob = avg_kras_allele_prob_weighted) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(barplot = purrr::map2(
        cancer, data, single_kras_allele_freq_barplot,
        y_max = max(barplot_df$avg_kras_allele_prob) * 1.1,
        with_errorbars = FALSE
    ))

predict_kras_allele_frequency_weighted_barplot2 <- patchwork::wrap_plots(
    predict_kras_allele_frequency_weighted_barplot2_df$barplot,
    nrow = 2
)

ggsave_wrapper(
    predict_kras_allele_frequency_weighted_barplot2,
    plot_path(GRAPHS_DIR,
              "predict_kras_allele_frequency_weighted_barplot2.svg"),
    "wide"
)

# The same bar plot as above, but separated by allele to compare across cancer.
predict_kras_allele_frequency_barplot3 <- barplot_df %>%
    ggplot(aes(
        x = cancer, y = avg_kras_allele_prob,
        color = cancer, fill = cancer
    )) +
    facet_wrap(~ kras_allele, scales = "free") +
    geom_col(alpha = 0.5) +
    geom_errorbar(
        aes(ymax = ci_upper_avg_kras_allele_prob,
            ymin = ci_lower_avg_kras_allele_prob),
        width = 0.2
    ) +
    scale_fill_manual(values = cancer_palette) +
    scale_color_manual(values = cancer_palette) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.02))) +
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
    predict_kras_allele_frequency_barplot3,
    plot_path(GRAPHS_DIR, "predict_kras_allele_frequency_barplot3.svg"),
    size = "large"
)



#### ---- Scatter plots ---- ####


# Make a single scatter plot for an individual cancer.
# This was made to be called from `obs_v_pred_scatter_plot` for each cancer.
individual_obs_v_pred_scatter_plot <- function(plot_df, cancer,
                                               stats = FALSE,
                                               p_value_limits = NULL,
                                               with_errorbars = FALSE) {

    if (with_errorbars) {
        axis_lim <- max(c(plot_df$real_kras_allele_frequency,
                          plot_df$ci_upper_avg_kras_allele_prob))
        plot_df %<>%
            mutate(ci_lower_avg_kras_allele_prob = purrr::map_dbl(
                ci_lower_avg_kras_allele_prob, ~ max(.x, 0)
            ))
    } else {
        axis_lim <- max(c(plot_df$real_kras_allele_frequency,
                          plot_df$avg_kras_allele_prob))
    }

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
                limits = p_value_limits,
                guide = guide_legend(
                    title.position = "left",
                    keywidth = unit(4, "mm"),
                    keyheight = unit(4, "mm"),
                    label.position = "right",
                    ncol = 1,
                    order = 1
            ))

    } else {
        p <- p + geom_point(aes(color = kras_allele),
                            size = 1, alpha = 0.7)
    }

    if (with_errorbars) {
        p <- p + geom_linerange(
                aes(ymin = ci_lower_avg_kras_allele_prob,
                    ymax = ci_upper_avg_kras_allele_prob,
                    color = kras_allele),
                alpha = 0.3
            )
    }

    p <- p +
        geom_abline(intercept = 0, slope = 1) +
        ggrepel::geom_text_repel(
            aes(label = kras_allele),
            color = "grey25",
            family = "Arial",
            size = 2,
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
            expand = expansion(mult = c(0, 0.02))
        ) +
        scale_y_continuous(
            limits = c(0, axis_lim),
            expand = expansion(mult = c(0, 0.02))
        ) +
        theme_bw(base_size = 8, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            legend.spacing = unit(c(10, 10), "mm"),
            legend.title = element_markdown(angle = 90,
                                            vjust = 0.5,
                                            hjust = 0.5,
                                            size = 6,
                                            family = "Arial"),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2, "mm"),
            legend.direction = "vertical",
            strip.background = element_blank()
        ) +
        labs(
            title = cancer,
            color = "KRAS allele",
            size = "-*log*<sub>10</sub>( p-value )",
            x = "observed frequency",
            y = "predicted frequency"
        )
    return(p)
}


# Scatter plot of observed vs. predicted KRAS allele frequency.
obs_v_pred_scatter_plot <- function(real_tib, pred_tib,
                                    save_template,
                                    stats_tib = NULL,
                                    p_value_cut = 0.05,
                                    with_errorbars = FALSE,
                                    save_size = "small") {
    pdata <- left_join(real_tib, pred_tib,
                        by = c("cancer", "kras_allele")) %>%
        filter(kras_allele %in% names(short_allele_pal)) %>%
        filter(!is.na(avg_kras_allele_prob)) %>%
        filter(cancer != "SKCM")

    p_val_range <- NULL
    if (!is.null(stats_tib)) {
        pdata <- stats_tib %>%
            select(cancer, kras_allele, p_value) %>%
            right_join(pdata, by = c("cancer", "kras_allele")) %>%
            mutate(log_p_value = -log10(p_value),
                   point_shape = ifelse(p_value < 0.05, 16, 17))

        p_val_range <- purrr::map_dbl(list(min, max),
                                      ~ .x(-log10(pdata$p_value)))
        p_val_range[[1]] <- max(floor(p_val_range[[1]] / 10) * 10, 0)
        p_val_range[[2]] <- ceiling(p_val_range[[2]] / 10) * 10
    }


    make_scatter_plot_and_save <- function(CANCER) {
        p <- pdata %>%
            filter(cancer == !!CANCER) %>%
            individual_obs_v_pred_scatter_plot(cancer = CANCER,
                                               stats = !is.null(stats_tib),
                                               p_value_limits = p_val_range,
                                               with_errorbars = with_errorbars)
        ggsave_wrapper(
            p,
            plot_path(GRAPHS_DIR, as.character(glue(save_template))),
            save_size
        )
        return(p)
    }

    all_plots <- purrr::map(
        sort(unique(pdata$cancer)), make_scatter_plot_and_save
    )

    invisible(all_plots)
}


# Scatter plot of all alleles observed vs. predicted KRAS allele frequency.
real_af <- cancer_muts_df %>%
    filter(ras == "KRAS") %>%
    mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
    remove_alleles_with_low_frequency() %>%
    select(cancer, tumor_sample_barcode, kras_allele) %>%
    group_by(cancer, tumor_sample_barcode, kras_allele) %>%
    slice(1) %>%
    ungroup() %>%
    unique() %>%
    count(cancer, kras_allele, name = "kras_allele_count") %>%
    group_by(cancer) %>%
    mutate(real_kras_allele_frequency = kras_allele_count / sum(kras_allele_count)) %>%
    ungroup()

predicted_af <- kras_hotspot_probability %>%
    filter(!is.na(kras_allele_prob)) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
    ) %>%
    ungroup() %>%
    left_join(kras_allele_freq_bootstrap_ci,
              by = c("cancer", "kras_allele")) %>%
    mutate(
        ci_upper_avg_kras_allele_prob = avg_kras_allele_prob + ci_upper,
        ci_lower_avg_kras_allele_prob = avg_kras_allele_prob - ci_lower
    )

obs_v_pred_scatter_plot(real_af, predicted_af,
                        save_template = "obs_pred_plot_{CANCER}.svg",
                        with_errorbars = TRUE)

# left_join(real_af, predicted_af,
#                     by = c("cancer", "kras_allele")) %>%
#     filter(kras_allele %in% names(short_allele_pal)) %>%
#     filter(!is.na(avg_kras_allele_prob)) %>%
#     filter(cancer != "SKCM")


# Scatter plot of G12 alleles observed vs. predicted KRAS allele frequency.
predicted_af_g12 <- predicted_kras_allele_frequency %>%
    filter(kras_codon == 12) %>%
    remove_alleles_with_low_frequency() %>%
    group_by(tumor_sample_barcode) %>%
    mutate(
        kras_allele_prob = tricontext_mut_count / tricontext_count,
        kras_allele_prob = kras_allele_prob / sum(kras_allele_prob)
    ) %>%
    filter(!is.na(kras_allele_prob)) %>%
    group_by(cancer, kras_allele) %>%
    summarise(
        avg_kras_allele_prob = mean(kras_allele_prob),
    ) %>%
    ungroup() %>%
    left_join(kras_allele_freq_bootstrap_ci,
              by = c("cancer", "kras_allele")) %>%
    mutate(
        ci_upper_avg_kras_allele_prob = avg_kras_allele_prob + ci_upper,
        ci_lower_avg_kras_allele_prob = avg_kras_allele_prob - ci_lower
    )

obs_v_pred_scatter_plot(real_af, predicted_af_g12,
                        save_template = "obs_pred_plot_g12_{CANCER}.svg",
                        with_errorbars = TRUE)




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
    b <- binom.test(x = df$kras_allele_count[[1]],
                    n = df$num_cancer_samples[[1]],
                    p = df$avg_kras_allele_prob[[1]])
    return(tidy(b))
}


kras_allele_freq_stats <- left_join(
        real_af, predicted_af, by = c("cancer", "kras_allele")
    ) %>%
    filter(!is.na(avg_kras_allele_prob)) %>%
    filter(kras_allele %in% names(short_allele_pal)) %>%
    filter(cancer != "SKCM") %>%
    left_join(cancer_sample_count_df, by = "cancer") %>%
    group_by(cancer, kras_allele) %>%
    nest() %>%
    ungroup() %>%
    mutate(binom_res = purrr::map(data, kras_allele_binomial_test)) %>%
    unnest(c(data, binom_res)) %>%
    janitor::clean_names()

kras_allele_freq_stats %>%
    select(-ci_obj) %T>%
    write_tsv(table_path(GRAPHS_DIR, "kras_allele_freq_stats.tsv")) %>%
    select(cancer, kras_allele,
           real_kras_allele_frequency, avg_kras_allele_prob,
           ci_lower_avg_kras_allele_prob, ci_upper_avg_kras_allele_prob,
           p_value, conf_low, conf_high, method, alternative) %>%
    knitr::kable(digits = 3)


kras_allele_freq_stats %>%
    filter(p_value >= 0.05) %>%
    select(cancer, kras_allele,
           real_kras_allele_frequency, avg_kras_allele_prob,
           p_value)

plots <- obs_v_pred_scatter_plot(
    real_af, predicted_af,
    stats_tib = kras_allele_freq_stats,
    with_errorbars = TRUE,
    save_template = "obs_pred_plot_stats_{CANCER}.svg"
)
saveFigRds(plots, "obs_pred_plot_stats")



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

kras_g12_freq_stats %>%
    select(-ci_obj) %>%
    write_tsv(table_path(GRAPHS_DIR, "kras_g12_freq_stats.tsv"))

kras_g12_freq_stats %>%
    filter(p_value >= 0.05) %>%
    select(cancer, kras_allele,
           real_kras_allele_frequency, avg_kras_allele_prob,
           p_value)

plots <- obs_v_pred_scatter_plot(real_af, predicted_af_g12,
                        stats_tib = kras_g12_freq_stats,
                        with_errorbars = TRUE,
                        save_template = "obs_pred_plot_g12_stats_{CANCER}.svg")
saveFigRds(plots, "obs_pred_plot_g12_stats")
