# Plot the distribution of mutational signatures by allele.

GRAPHS_DIR <- "50_20_mutsignatures-distributions"
TABLES_DIR <- "50_20_mutsignatures-distributions"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(TABLES_DIR)

alleles_frequency_per_cancer_df <- mutational_signatures_df %>%
    filter(!is_hypermutant) %>%
    mutate(ras_allele = str_remove_all(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, cancer, ras_allele) %>%
    unique() %>%
    count(cancer, ras_allele)

# Get the alleles that reach a minimum number `min_num` of samples.
alleles_to_plot <- function(cancer, min_num = 15) {
    alleles_frequency_per_cancer_df %>%
        filter(cancer == !!cancer & n >= min_num) %>%
        pull(ras_allele) %>%
        unique()
}



#### ---- Mutational signature levels for all samples ---- ####


# Pivot the long tibble to a wide data frame for h. clustering.
data_to_dataframe_for_clustering <- function(dat) {
    df <- dat %>%
        ungroup() %>%
        select(tumor_sample_barcode, description, contribution) %>%
        pivot_wider(names_from = description,
                    values_from = contribution) %>%
        as.data.frame() %>%
        column_to_rownames("tumor_sample_barcode")
    return(df)
}


# Perform hierarchical clustering on a wide data frame.
# Prepare `mutsig_noartifact_df` with `data_to_dataframe_for_clustering()`.
dist_hclust <- function(df) {
    hc <- df %>%
        dist(method = "manhattan") %>%
        hclust(method = "complete")
    return(hc)
}


# Perform scaling and hierarchical clustering on a wide data frame.
# Prepare `mutsig_noartifact_df` with `data_to_dataframe_for_clustering()`.
scale_dist_hclust <- function(df) {
    scaled_df <- as.data.frame(t(apply(df, 1, scale)))
    colnames(scaled_df) <- colnames(df)
    return(dist_hclust(scaled_df))
}


# Order the samples by hierarchical clustering.
# Input should have columns:
#     'tumor_sample_barcode', 'description', 'contribution'
# Output is the same as the input data frame with a new column 'sample_idx'.
order_samples_by_hclust <- function(dat) {
    hc <- data_to_dataframe_for_clustering(dat) %>%
        dist_hclust()

    sample_order_df <- tibble(
        tumor_sample_barcode = hc$labels,
        sample_idx = hc$order
    )

    mod_dat <- full_join(dat, sample_order_df, by = "tumor_sample_barcode")
    if(any(is.na(mod_dat))) {
        stop("Joining in `order_samples_by_hclust()` went wrong.")
    }
    return(mod_dat)
}


# Order the samples by MEMO sort.
# Input should have columns:
#     'tumor_sample_barcode', 'description', 'contribution'
# Output is the same as the input data frame with a new column 'sample_idx'.
order_samples_by_memosort <- function(dat) {

    if (!any(colnames(dat) == "signature_idx")) {
        stop("Sort signatures before calling `order_samples_by_memosort()`.")
    }

    sig_score_tib <- dat %>%
        select(description, signature_idx, contribution) %>%
        unique() %>%
        mutate(value = (2 ^ (signature_idx)) - 1) %>%
        select(-signature_idx, -contribution)


    sample_order_df <- dat %>%
        select(tumor_sample_barcode, description, contribution) %>%
        unique() %>%
        left_join(sig_score_tib, by = "description") %>%
        mutate(value = value * contribution) %>%
        group_by(tumor_sample_barcode) %>%
        summarise(
            score = sum(value)
        ) %>%
        ungroup() %>%
        arrange(score) %>%
        mutate(sample_idx = seq(1:n())) %>%
        select(-score)

    mod_dat <- full_join(dat, sample_order_df, by = "tumor_sample_barcode")
    if(any(is.na(mod_dat))) {
        stop("Joining in `order_samples_by_hclust()` went wrong.")
    }
    return(mod_dat)
}


# Order the signatures by hierarchical clustering.
# Input should have columns:
#     'tumor_sample_barcode', 'description', 'contribution'
# Output is the same as the input data frame with a new column 'signature_idx'.
order_signatures_by_hclust <- function(dat) {
    hc <- data_to_dataframe_for_clustering(dat) %>%
        t() %>%
        dist_hclust()

    sig_order_df <- tibble(
        description = hc$labels,
        signature_idx = hc$order
    )

    mod_dat <- full_join(dat, sig_order_df, by = "description")
    if(any(is.na(mod_dat))) {
        stop("Joining in `order_signatures_by_hclust()` went wrong.")
    }
    return(mod_dat)
}


# Order the signatures by mean value.
# Input should have columns:
#     'tumor_sample_barcode', 'description', 'contribution'
# Output is the same as the input data frame with a new column 'signature_idx'.
order_signatures_by_avg <- function(dat) {
    sig_order_df <- dat %>%
        group_by(description) %>%
        summarise(avg_val = mean(description)) %>%
        ungroup() %>%
        arrange(avg_val) %>%
        mutate(signature_idx = seq(1, n())) %>%
        select(description, signature_idx)

    mod_dat <- full_join(dat, sig_order_df, by = "description")
    if(any(is.na(mod_dat))) {
        stop("Joining in `order_signatures_by_avg()` went wrong.")
    }
    return(mod_dat)
}


# Prepare the data for bar-plots with `barplot_distribution_per_sample()`.
# Option to order the signatures.
prep_barplot_distribution_per_sample <- function(data,
                                                 order_signatures = FALSE) {
    plot_data <- data
    plot_data$tumor_sample_barcode <- fct_reorder(
        plot_data$tumor_sample_barcode,
        plot_data$sample_idx
    )

    if (order_signatures) {
        plot_data$description <- fct_reorder(
            plot_data$description,
            plot_data$signature_idx
        )
    }

    return(plot_data)
}


# Stacked bar-plot of the mutational signatures for each tumor sample.
barplot_distribution_per_sample <- function(data, title) {
    p <- data %>%
        ggplot(aes(x = tumor_sample_barcode, y = contribution)) +
        geom_col(
            aes(fill = description),
            position = "fill",
            width = 1.0
        ) +
        scale_fill_manual(
            values = mutsig_descrpt_pal,
            guide = guide_legend(
                nrow = 1,
                label.position = "top",
                label.hjust = 0.5,
                label.vjust = -2
            )
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_bw(
            base_size = 7,
            base_family = "Arial"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.key.size = unit(4, "mm")
        ) +
        labs(
            title = title,
            x = "tumor samples",
            y = "level",
            fill = "signature"
        )
    return(p)
}


# Save a plot to `GRAPHS_DIR`.
save_barplot_distribution_per_sample_plots <- function(cancer, gg_obj, ...) {
    ggsave_wrapper(
        gg_obj,
        plot_path(GRAPHS_DIR, glue("signature-level-per-sample_{cancer}.jpeg")),
        width = 12, height = 2
    )
    return(gg_obj)
}



# Stacked bar-plot of mutational signatures per tumor sample barcode.
# One plot per cancer.
mutsig_per_sample_plots <- mutsig_noartifact_df %>%
    group_by(tumor_sample_barcode) %>%
    filter(!all(contribution == 0)) %>%
    group_by(cancer, description) %>%
    filter(!all(contribution == 0)) %>%
    group_by(tumor_sample_barcode, cancer, description) %>%
    summarise(contribution = sum(contribution)) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        data = purrr::map(data, order_signatures_by_avg),
        data = purrr::map(data, order_samples_by_memosort),
        data = purrr::map(data,
                          prep_barplot_distribution_per_sample,
                          order_signatures = FALSE),
        gg_obj = purrr::map2(data, cancer, barplot_distribution_per_sample)
    ) %>%
    pwalk(save_barplot_distribution_per_sample_plots)

# Save faceted plot as proto for Supp. Figure 2.
saveFigRds(mutsig_per_sample_plots, "mutsig_per_sample_plots")



# Stacked bar-plot of mutational signatures per tumor sample barcode.
# Facet by cancer.
mutsig_per_sample_plot <- mutsig_noartifact_df %>%
    filter(!all_zeros) %>%
    group_by(tumor_sample_barcode, cancer, description) %>%
    summarise(contribution = sum(contribution)) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        data = purrr::map(data, order_signatures_by_avg),
        data = purrr::map(data, order_samples_by_memosort),
        data = purrr::map(data,
                          prep_barplot_distribution_per_sample,
                          order_signatures = FALSE)
    ) %>%
    unnest(data) %>%
    barplot_distribution_per_sample(title = NULL) +
    facet_wrap(~cancer, scales = "free", ncol = 1)

# Save faceted plot as JPEG.
ggsave_wrapper(
    mutsig_per_sample_plot,
    plot_path(GRAPHS_DIR, "signature-level-per-sample.jpeg"),
    "large"
)
# Save faceted plot as SVG
ggsave_wrapper(
    mutsig_per_sample_plot,
    plot_path(GRAPHS_DIR, "signature-level-per-sample.svg"),
    "large"
)
# Save faceted plot as proto for Supp. Figure 2.
saveFigRds(mutsig_per_sample_plot, "signature-level-per-sample")



#### ---- Boxplots of mutational signature levels ---- ####

# Make boxplots of the signature levels in all samples.
signature_distribution_boxplots <- function(tib) {
    p <- tib %>%
        mutate(
            description = factor(description,
                                 levels = names(mutsig_descrpt_pal))
        ) %>%
        ggplot(aes(x = description, y = contribution)) +
        facet_wrap(~ cancer, scales = "free", ncol = 1) +
        geom_boxplot(
            aes(fill = description),
            outlier.shape = NA
        ) +
        scale_fill_manual(
            values = mutsig_descrpt_pal
        ) +
        theme_bw(
            base_size = 7,
            base_family = "Arial"
        ) +
        labs(
            x = "mutational signature",
            y = "level"
        )
    return(p)
}

# Boxplots of mutational signatures in each sample (values of zero ignored).
sig_boxes_with0s <- mutsig_noartifact_df %>%
    group_by(tumor_sample_barcode, cancer, description) %>%
    summarise(contribution = sum(contribution)) %>%
    ungroup() %>%
    filter(contribution > 0) %>%
    signature_distribution_boxplots()
ggsave_wrapper(
    sig_boxes_with0s,
    plot_path(GRAPHS_DIR, "signature-level-boxplots_with0.svg"),
    "tall"
)
saveFigRds(sig_boxes_with0s, "signature-level-boxplots_with0")

# Boxplots of mutational signatures in each sample (values of zero maintained).
sig_boxes <- mutsig_noartifact_df %>%
    group_by(tumor_sample_barcode, cancer, description) %>%
    summarise(contribution = sum(contribution)) %>%
    ungroup() %>%
    signature_distribution_boxplots()
ggsave_wrapper(
    sig_boxes,
    plot_path(GRAPHS_DIR, "signature-level-boxplots.svg"),
    "tall"
)
saveFigRds(sig_boxes, "signature-level-boxplots")



#### ---- Mutational signatures per KRAS allele ---- ####

# Combine the plots into a grid using 'patchwork'.
distribution_plots <- mutational_signatures_df %>%
    filter(description != "artifact" & !is_hypermutant) %>%
    mutate(
        signature = factor(signature, levels = names(mutsig_pal)),
        description = factor(description, levels = names(mutsig_descrpt_pal)),
        ras_allele = str_remove_all(ras_allele, "KRAS_"),
        ras_allele = factor_alleles(ras_allele)
    ) %>%
    group_by(cancer) %>%
    filter(ras_allele %in% alleles_to_plot(cancer = unique(cancer))) %>%
    group_by(cancer, ras_allele, description) %>%
    summarise(contribution = sum(contribution)) %>%
    ungroup() %>%
    ggplot(aes(x = ras_allele, y = contribution)) +
    facet_wrap(~ cancer, scales = "free", nrow = 2) +
    geom_col(aes(fill = description), position = "fill") +
    scale_fill_manual(
        values = mutsig_descrpt_pal,
        guide = guide_legend(
            title = "signature",
            nrow = 1,
            title.position = "left",
            title.vjust = 0.5,
            label.position = "top",
            label.hjust = 0.5,
            label.vjust = -6
        )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw(
        base_size = 7,
        base_family = "Arial"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(4, "mm"),
        legend.spacing.y = unit(0, "mm"),
        strip.background = element_blank()
    ) +
    labs(
        y = "average level"
    )

ggsave_wrapper(
    distribution_plots,
    plot_path(GRAPHS_DIR, "mutsig-dist_combined.svg"),
    'wide'
)

saveFigRds(distribution_plots,
           "mutational-signatures-distribution-by-allele")



#### ---- Clock vs. non-clock ---- ####

clock_sigs <- c("1", "5")
clock_pal <- c(
    "clock" = "grey25",
    "non-clock" = "grey80"
)

# A theme for the clock vs. non-clock plots (below)
clock_plot_theme <- function() {
    theme_bw(
        base_size = 7,
        base_family = "Arial"
    ) %+replace%
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none",
            legend.title = element_blank(),
            strip.background = element_blank()
        )
}

# Prepare the signatures to plot clock vs. non-clock.
prepare_clock_signature_data <- function(data) {
    plot_data <- data %>%
        filter(!is_hypermutant) %>%
        mutate(clock_signature = ifelse(
            signature %in% !!clock_sigs, "clock", "non-clock"
        )) %>%
        group_by(cancer, tumor_sample_barcode, clock_signature) %>%
        summarise(contribution = sum(contribution)) %>%
        ungroup()
    return(plot_data)
}

# Distribution of levels of clock vs. non-clock signatures.
# A violin with overlaid boxplots.
clock_violin_box <- mutsig_noartifact_df %>%
    prepare_clock_signature_data() %>%
    ggplot(aes(x = clock_signature, y = contribution)) +
    geom_violin(
        aes(fill = clock_signature),
        alpha = 0.7
    ) +
    geom_boxplot(
        fill = "white",
        width = 0.15,
        outlier.shape = NA,
        notch = FALSE
    ) +
    scale_fill_manual(
        values = clock_pal
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    clock_plot_theme() +
    facet_wrap(~ cancer, scales = "free", nrow = 2) +
    labs(
        y = "total level"
    )

ggsave_wrapper(
    clock_violin_box,
    plot_path(GRAPHS_DIR, "clock-signatures_violin-box.svg"),
    width = 5, height = 4
)

saveFigRds(clock_violin_box,
           "clock-signatures_violin-box")

# Distribution of levels of clock vs. non-clock signatures.
# horizontal density plot.
clock_ridges <- mutsig_noartifact_df %>%
    prepare_clock_signature_data() %>%
    ggplot(aes(y = clock_signature, x = contribution)) +
    ggridges::geom_density_ridges(
        aes(fill = clock_signature),
        alpha = 0.7
    ) +
    scale_fill_manual(
        values = clock_pal
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    clock_plot_theme() +
    facet_wrap(~ cancer, scales = "free", nrow = 2) +
    labs(
        y = "total level"
    )

ggsave_wrapper(
    clock_ridges,
    plot_path(GRAPHS_DIR, "clock-signatures_ridges.svg"),
    width = 5, height = 4
)




#### ---- Plot Clock vs. Smoke contribution in KRAS mutants ---- ####

clock_vs_smoke_plot <- mutsig_noartifact_df %>%
    filter(
        cancer == "LUAD" &
        signature %in% c(clock_sigs, "4"),
        ras_allele != "WT"
    ) %>%
    mutate(signature = paste0("Sig_", signature)) %>%
    select(tumor_sample_barcode, ras_allele, signature, contribution) %>%
    pivot_wider(
        names_from = signature,
        values_from = contribution
    ) %>%
    mutate(
        clock = Sig_1 + Sig_5,
        smoke = Sig_4,
        allele = str_remove(ras_allele, "KRAS_"),
        allele = ifelse(
            allele %in% names(short_allele_pal), allele, "Other"
        ),
        allele = factor_alleles(allele),
        alpha = ifelse(allele == "G12C", 0.8, 0.6),
        size = ifelse(allele == "G12C", 1.1, 0.7)
    ) %>%
    ggplot(aes(x = clock, y = smoke)) +
    geom_point(
        aes(color = allele, alpha = alpha, size = size)
    ) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
    scale_color_manual(values = short_allele_pal) +
    scale_alpha_identity() +
    scale_size_identity() +
    theme_bw(
        base_size = 8,
        base_family = "Arial"
    ) +
    theme(
        legend.title = element_blank(),
    ) +
    labs(
        title = "Clock vs. smoke in LUAD of KRAS mutants.",
        subtitle = "There is a focus on the smoking-related allele, G12C."
    )
ggsave_wrapper(
    clock_vs_smoke_plot,
    plot_path(GRAPHS_DIR, "clock_vs_smoke.svg"),
    "small"
)


#### ---- Tabkes of signature levels ---- ####

# Table of all mutational signatures for each sample.
mutational_signatures_df %>%
    mutate(allele = str_remove_all(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, signature, contribution, dataset, target,
           cancer, allele, is_hypermutant) %>%
    pivot_wider(names_from = signature,
                names_prefix = "sig_",
                values_from = contribution) %>%
    write_tsv(table_path(TABLES_DIR, "mutation-signature-contribution_all.tsv"))


# Table of all mutational signatures for each sample without "Artifact" sig.
mutsig_noartifact_df %>%
    mutate(allele = str_remove_all(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, signature, contribution, dataset, target,
           cancer, allele, is_hypermutant) %>%
    pivot_wider(names_from = signature,
                names_prefix = "sig_",
                values_from = contribution) %>%
    write_tsv(table_path(TABLES_DIR, "mutation-signature-contribution_no-artifact.tsv"))



# Smoking vs. clock cumulative levels in each sample
clock_vs_smoke_df <- mutsig_noartifact_df %>%
    filter(cancer == "LUAD" & signature %in% c(clock_sigs, "4")) %>%
    mutate(allele = str_remove(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, allele, signature, contribution) %>%
    pivot_wider(names_from = signature,
                names_prefix = "sig_",
                values_from = contribution) %>%
    mutate(smoke = sig_4,
           clock = sig_1 + sig_5)

# Table of levels for each sample.
clock_vs_smoke_df %>%
    arrange(allele, -smoke, clock) %>%
    write_tsv(table_path(TABLES_DIR, "smoke-vs-clock_LUAD_all.tsv"))


# Summarized values for each KRAS allele.
clock_vs_smoke_df %>%
    group_by(allele) %>%
    summarise(
        num_cancer_samples = n_distinct(tumor_sample_barcode),
        avg_smoke = mean(smoke),
        avg_clock = mean(clock),
        sd_smoke = sd(smoke),
        sd_clock = sd(clock)
    ) %>%
    ungroup() %>%
    arrange(-num_cancer_samples) %>%
    write_tsv(table_path(TABLES_DIR, "smoke-vs-clock_LUAD_summary.tsv"))
