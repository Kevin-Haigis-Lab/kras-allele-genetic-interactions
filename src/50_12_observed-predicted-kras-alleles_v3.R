# Observed vs. predicted KRAS allele frequency using the trinucleotide context.
# This is a new version for "src/50_10_observed-predicted-kras-alleles.R" because
# I am not confident that it is bug-free.
#
# Also, GM and I discussed the statsitical testing and decided to calculate and
# R-squared value and use the Chi squared test instead of a binomial.

set.seed(0)

GRAPHS_DIR <- "50_12_observed-predicted-kras-alleles_v3"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)

mut_sig_descriptions <- mutational_signatures_df %>%
    distinct(signature, description)

artifact_signatures <- mut_sig_descriptions %>%
    filter(description == "artifact") %>%
    mutate(signature = paste0("sig", signature)) %>%
    pull(signature)

msg_sigs <- paste0(artifact_signatures, collapse = ", ")
message(glue("The following are artifact signatures: {msg_sigs}"))


oncogenic_alleles <- unique(kras_trinucleotide_contexts$kras_allele)
message("The following are the oncogenic KRAS alleles:")
print(oncogenic_alleles)


#### ---- KRAS allele frequency ---- ####

# Calculate the frequency of the KRAS mutant alleles.
# This only includes the oncogenic KRAS alleles.
calc_frequency_of_oncogenic_kras <- function(alleles) {
    allele_tbl <- table(alleles)
    allele_tbl <- allele_tbl[names(allele_tbl) %in% oncogenic_alleles]
    allele_tbl <- as_tibble(allele_tbl) %>%
        mutate(n = n / sum(n)) %>%
        set_names(c("kras_allele", "frequency"))

    missing_idx <- !(oncogenic_alleles %in% allele_tbl$kras_allele)
    missing_alleles <- oncogenic_alleles[missing_idx]

    if (length(missing_alleles) > 0) {
        missing_tbl <- tibble(kras_allele = missing_alleles,
                              frequency = 0.0)
        allele_tbl <- bind_rows(allele_tbl, missing_tbl)
    }
    return(allele_tbl)
}


real_kras_allele_freq <- cancer_full_coding_muts_df %>%
    filter(cancer != "SKCM") %>%
    distinct(cancer, tumor_sample_barcode, ras_allele) %>%
    mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
    group_by(cancer) %>%
    summarise(
        kras_freq = list(calc_frequency_of_oncogenic_kras(kras_allele))
    ) %>%
    unnest(kras_freq) %>%
    dplyr::rename(real_allele_freq = frequency)


#### ---- Samples with too few mutations ---- ####

MIN_NUM_MUTATION_PER_SAMPLE <- 20

REMOVE_TSB <- trinucleotide_mutations_df %>%
    filter(cancer != "SKCM") %>%
    count(tumor_sample_barcode) %>%
    filter(n < MIN_NUM_MUTATION_PER_SAMPLE) %>%
    pull(tumor_sample_barcode)

message(glue(
    "Removing {length(REMOVE_TSB)} samples because they have too few mutations"
))


#### ---- Samples with all artifact mutational signatures ---- ####

MAXIMUM_ARTIFACT_SIGNATURE <- 0.50
artifact_contribution <- mutational_signatures_df %>%
    filter(cancer != "SKCM") %>%
    mutate(signature = paste0("sig", signature)) %>%
    filter(!tumor_sample_barcode %in% REMOVE_TSB &
           signature %in% artifact_signatures) %>%
    group_by(tumor_sample_barcode) %>%
    summarise(contribution = sum(contribution)) %>%
    ungroup()

artifact_contribution_plt <- artifact_contribution %>%
    filter(contribution > 0) %>%
    ggplot(aes(x = contribution)) +
    geom_histogram(bins = 100) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "contribution of Artifact",
         y = "count")
ggsave_wrapper(
    artifact_contribution_plt,
    plot_path(GRAPHS_DIR, "artifact_contribution_plt.svg"),
    "small"
)

ARTIFACT_TSB <- artifact_contribution %>%
    filter(contribution > !!MAXIMUM_ARTIFACT_SIGNATURE) %>%
    pull(tumor_sample_barcode) %>%
    unique()

message(glue(
    "Removing {length(ARTIFACT_TSB)} samples because they have too much Artifact signature."
))

REMOVE_TSB <- c(REMOVE_TSB, ARTIFACT_TSB)

# Remaining number of samples per cancer.
trinucleotide_mutations_df %>%
    filter(cancer != "SKCM") %>%
    filter(!tumor_sample_barcode %in% REMOVE_TSB) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    select(cancer, tumor_sample_barcode) %>%
    count(cancer, name = "num_tumor_samples") %>%
    knitr::kable()

#> |cancer | num_tumor_samples|
#> |:------|-----------------:|
#> |COAD   |              1500|
#> |LUAD   |               852|
#> |MM     |              1185|
#> |PAAD   |              1192|


tricontext_genome_counts <- tricontext_counts_df %>%
    set_names(c("context", "exome", "genome")) %>%
    pivot_longer(-context, names_to = "target", values_to = "genome_count")


MINIMUM_ALLELE_FREQ <- 0.03
TOP_ALLELES_PER_CANCER <- 5
alleles_for_each_cancer <- real_kras_allele_freq %>%
    filter(cancer != "SKCM") %>%
    arrange(cancer, desc(real_allele_freq)) %>%
    group_by(cancer) %>%
    mutate(cancer_rank = 1:n()) %>%
    ungroup() %>%
    filter(real_allele_freq >= MINIMUM_ALLELE_FREQ |
           cancer_rank <= TOP_ALLELES_PER_CANCER) %>%
    select(cancer, kras_allele)

alleles_for_each_cancer %>%
    group_by(cancer) %>%
    summarise(alleles = paste0(kras_allele, collapse = ", "))


kras_allele_predictions <- trinucleotide_mutations_df %>%
    filter(cancer != "SKCM" & !(tumor_sample_barcode %in% REMOVE_TSB)) %>%
    count(cancer, tumor_sample_barcode, target, context, tricontext,
          name = "tumor_count") %>%
    left_join(tricontext_genome_counts, by = c("context", "target")) %>%
    mutate(tumor_count_norm = tumor_count / genome_count) %>%
    inner_join(kras_trinucleotide_contexts,
               by = c("context", "tricontext")) %>%
    right_join(alleles_for_each_cancer, by = c("cancer", "kras_allele")) %>%
    group_by(cancer, tumor_sample_barcode, target, kras_allele) %>%
    summarise(allele_prob = sum(tumor_count_norm)) %>%
    group_by(cancer, tumor_sample_barcode) %>%
    mutate(allele_prob = allele_prob / sum(allele_prob)) %>%
    ungroup() %>%
    left_join(real_kras_allele_freq, by = c("cancer", "kras_allele"))

# KRAS allele probabilities per tumor sample.
kras_allele_predictions



#### ---- Statistics: boostrapping 95% CI ---- ####

# Calculate the expect frequency for the data of a single cancer.
calc_expected_frequency <- function(df) {
    df %>%
        group_by(kras_allele) %>%
        summarise(expected_allele_frequency = sum(allele_prob),
                  observed_allele_frequency = unique(real_allele_freq)) %>%
        ungroup() %>%
        mutate(
            expected_allele_frequency = expected_allele_frequency /
                                        sum(expected_allele_frequency),
            observed_allele_frequency = observed_allele_frequency /
                                        sum(observed_allele_frequency)
        ) %>%
        ungroup()
}


# Bootstrap CIs for the expected frequencies.
calc_expected_frequency_boot <- function(data,
                                         index = NULL,
                                         all_alleles = NULL) {
    exp_freqs <- data[index, ] %>%
        unnest(data) %>%
        calc_expected_frequency() %>%
        select(-observed_allele_frequency)

    allele_results <- tibble(kras_allele = all_alleles) %>%
        left_join(exp_freqs, by = "kras_allele") %>%
        mutate(expected_allele_frequency = ifelse(
            is.na(expected_allele_frequency), 0, expected_allele_frequency
        ))

    return(deframe(allele_results))
}


# Real values
cancer_expect_frequencies <- kras_allele_predictions %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(data = map(data, calc_expected_frequency)) %>%
    unnest(data)


# Bootstrapping
boot_cancer_expect_frequncies <- function(cancer, df, R = 1e3) {
    nested_df <- df %>% group_by(tumor_sample_barcode) %>% nest()
    boot_res <- boot::boot(nested_df,
               calc_expected_frequency_boot,
               R = R,
               all_alleles = unique(df$kras_allele))
    return(list(boot_obj = boot_res,
                alleles = unique(df$kras_allele)))
}


# Extract the CIs from a boot object at an index.
get_conf_intervals <- function(boot_obj, conf, index) {
    boot::boot.ci(boot_obj,
                  index = index,
                  conf = conf,
                  type = "basic")$basic[, c(4, 5)]
}

# Extract 95% CI for each index from bootstrap results.
extract_boot_results <- function(boot_res, conf = 0.95) {
    f <- function(x, i) {
        ci <- get_conf_intervals(boot_res$boot_obj, conf = conf, index = i)
        tibble(kras_allele = x, lower_ci = ci[[1]], upper_ci = ci[[2]])
    }
    imap(boot_res$alleles, f) %>%
        bind_rows()
}


# Results from bootstrapping the samples used for the calculation of the
# predicted KRAS alleles.
cache("kras_allele_predictions_boot_results",
      depends = c("kras_allele_predictions"),
{
    set.seed(0)

    kras_allele_predictions_boot_results <- kras_allele_predictions %>%
        group_by(cancer) %>%
        nest() %>%
        mutate(
            boot_res = map2(cancer, data, boot_cancer_expect_frequncies,
                            R = 1e3),
            boot_ci = map(boot_res, extract_boot_results)
        )
})

cancer_expect_frequencies <- kras_allele_predictions_boot_results %>%
    select(cancer, boot_ci) %>%
    unnest(boot_ci) %>%
    right_join(cancer_expect_frequencies, by = c("cancer", "kras_allele"))





#### ---- Statistics: R-squared ---- ####

calc_obs_pred_rsquared <- function(x, y) {
    dat <- tibble(x = x, y = y)
    lm(y ~ x, data = dat) %>%
        broom::glance() %>%
        janitor::clean_names()
}


obs_pred_rsquareds <- cancer_expect_frequencies %>%
    group_by(cancer) %>%
    summarise(model_fit = list(
        calc_obs_pred_rsquared(observed_allele_frequency,
                               expected_allele_frequency)
    )) %>%
    unnest(model_fit) %>%
    ungroup() %>%
    select(cancer, r_squared, adj_r_squared, p_value) %>%
    rename(model_p_value = p_value)


#### ---- Statistics: Chi-Squared ---- ####


# Chi-squared test to nest null hypothesis that observed and predicted
# frequency of an allele are the same.
allele_pred_obs_chisquared <- function(num_mut, num_tot, pred_freq) {
    pred_alleles <- num_tot * pred_freq
    mat <- matrix(
        c(num_mut, num_tot - num_mut,
          pred_alleles, num_tot - pred_alleles),
        nrow = 2,
        byrow = TRUE
    )
    janitor::clean_names(broom::glance(chisq.test(mat)))
}


obs_pred_chisquared_res <- cancer_full_coding_muts_df %>%
    mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
    right_join(alleles_for_each_cancer, by = c("cancer", "kras_allele")) %>%
    group_by(cancer, kras_allele) %>%
    summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer) %>%
    mutate(num_cancer_samples = sum(num_allele_samples)) %>%
    ungroup() %>%
    left_join(cancer_expect_frequencies, by = c("cancer", "kras_allele")) %>%
    filter(!is.na(expected_allele_frequency)) %>%
    group_by(cancer, kras_allele) %>%
    mutate(
        chi_squared_test = list(
            allele_pred_obs_chisquared(num_allele_samples,
                                       num_cancer_samples,
                                       expected_allele_frequency)
        )
    ) %>%
    ungroup() %>%
    select(cancer, kras_allele, chi_squared_test) %>%
    unnest(chi_squared_test)

obs_pred_chisquared_res %>%
    filter(p_value > 0.05) %>%
    select(cancer, kras_allele, p_value) %>%
    knitr::kable(digits = 3)
#> |cancer |kras_allele | p_value|
#> |:------|:-----------|-------:|
#> |COAD   |G12A        |   0.637|
#> |LUAD   |G12V        |   0.896|
#> |MM     |G12A        |   0.863|
#> |MM     |G12D        |   0.424|
#> |MM     |G12R        |   0.420|
#> |MM     |G12V        |   0.689|
#> |MM     |G13D        |   0.238|
#> |MM     |Q61L        |   0.940|
#> |MM     |Q61R        |   0.519|


#### ---- Plot: Distribution of probabilities ---- ####

return_one_followed_by_NA <- function(x) {
    c(x[[1]], rep(NA, length(x) - 1))
}


per_sample_allele_probability_plot <- kras_allele_predictions %>%
    group_by(cancer, kras_allele) %>%
    mutate(real_allele_freq = return_one_followed_by_NA(real_allele_freq)) %>%
    ungroup() %>%
    ggplot(aes(x = allele_prob)) +
    facet_grid(cancer ~ kras_allele, scales = "free_y") +
    geom_density(aes(color = kras_allele), fill = NA) +
    geom_vline(aes(xintercept = real_allele_freq, color = kras_allele),
               lty = 2, alpha = 0.6) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, length.out = 5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_color_manual(values = short_allele_pal) +
    theme_bw(base_size = 7, base_family = "arial") +
    theme(
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 4.5, vjust = 1),
        axis.text.y = element_text(size = 4.5, hjust = 1),
        axis.ticks = element_blank()
    ) +
    labs(
        x = "probability per sample",
        y = "density",
        title = "Per sample probability of each KRAS allele",
        subtitle = "Separated by cancer; color relates to KRAS allele; vertical lines are the actual frequency"
    )
ggsave_wrapper(
    per_sample_allele_probability_plot,
    plot_path(GRAPHS_DIR, "per_sample_allele_probability_plot.svg"),
    "wide"
)


#### ---- Plot: Predicted vs. Observed ---- ####

plot_kras_allele_predictions <- function(cancer, data,
                                         p_val_cut = 0.05,
                                         use_cancer_color = FALSE) {

    pval_labels <- c(glue("p < {p_val_cut}"),
                     glue("p ≥ {p_val_cut}"))

    mod_data <- data %>%
        mutate(
            lower_ci = scales::squish(lower_ci, range = c(0, 1)),
            upper_ci = scales::squish(upper_ci, range = c(0, 1)),
            is_significant = ifelse(
                p_value < p_val_cut, pval_labels[[1]], pval_labels[[2]]
            ),
            is_significant = factor(is_significant, levels = pval_labels)
        )

    max_val <- max(c(mod_data$observed_allele_frequency,
                     mod_data$upper_ci))

    point_color <- ifelse(use_cancer_color,
                          cancer_palette[[cancer]],
                          "grey20")

    p <- mod_data %>%
        ggplot(aes(x = observed_allele_frequency,
                   y = expected_allele_frequency)) +
        geom_abline(lty = 2, size = 0.6, color = "grey60") +
        geom_linerange(aes(ymin = lower_ci, ymax = upper_ci),
                       color = "grey35",
                       size = 0.6) +
        geom_point(aes(shape = is_significant),
                   size = 1.3,
                   color = point_color) +
        ggrepel::geom_text_repel(aes(label = kras_allele),
                                 size = 2.2,
                                 family = "Arial",
                                 seed = 0,
                                 force = 1,
                                 box.padding = unit(2, "mm"),
                                 segment.alpha = 0.5,
                                 segment.color = "grey30",
                                 segment.size = 0.3,
                                 min.segment.length = unit(5, "mm")) +
        scale_x_continuous(limits = c(0, max_val),
                           expand = expansion(mult = c(0, 0.03))) +
        scale_y_continuous(limits = c(0, max_val),
                           expand = expansion(mult = c(0, 0.03)),
                           labels = function(x) {ifelse(x == 0, "", x)}) +
        scale_shape_manual(values = c(17, 16),
                           drop = FALSE,
                           guide = guide_legend(title.position = "top",
                                                label.position = "top")) +
        coord_fixed() +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            strip.background = element_blank(),
            legend.title = element_markdown(hjust = 0.5, vjust = 0.5)
        ) +
        labs(x = "observed",
             y = "predicted",
             shape = "χ<sup>2</sup> test",
             title = cancer)

    if (is.null(mod_data$adj_r_squared)) { return(p) }

    r_sq <- round(mod_data$adj_r_squared[[1]], 3)
    r_sq <- str_pad(r_sq, width = 5, side = "right", pad = "0")
    mdl_p_val <- round(mod_data$model_p_value[[1]], 2)
    mdl_p_val <- str_pad(mdl_p_val, width = 4, side = "right", pad = "0")
    lbl <- glue("adj. R<sup>2</sup> = {r_sq}<br>p-val = {mdl_p_val}")


    r_sq_xpos <- max_val * 0.05
    r_sq_ypos <- max_val - (max_val * 0.01)

    p +
        geom_richtext(
            label = lbl, x = r_sq_xpos, y = r_sq_ypos,
            hjust = 0, vjust = 1, family = "Arial", size = 2.8,
            fill = NA, label.color = NA,
            label.padding = grid::unit(rep(0, 4), "pt")
        )
}


save_kras_allele_predictions <- function(cancer, plt, gl_template) {
    ggsave_wrapper(
        plt,
        plot_path(GRAPHS_DIR, as.character(glue(gl_template))),
        size = "small"
    )

    saveFigRds(plt, as.character(glue(gl_template)))

    invisible(plt)
}


predicted_allele_frequency_scatter <- cancer_expect_frequencies %>%
    left_join(obs_pred_chisquared_res, by = c("cancer", "kras_allele")) %>%
    left_join(obs_pred_rsquareds, by = c("cancer")) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        plt = map2(cancer, data, plot_kras_allele_predictions,
                   use_cancer_color = TRUE),
        plt = map2(cancer, plt, save_kras_allele_predictions,
                   gl_template = "{cancer}_predict-allele-freq_scatter.svg")
    )
