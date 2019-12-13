# Plot the distribution of mutational signatures by allele.

GRAPHS_DIR <- "50_20_mutsignatures-distributions"
reset_graph_directory(GRAPHS_DIR)

alleles_frequency_per_cancer_df <- mutational_signatures_df %>%
    filter(!is_hypermutant) %>%
    mutate(ras_allele = str_remove_all(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, cancer, ras_allele) %>%
    unique() %>%
    count(cancer, ras_allele)

alleles_to_plot <- function(cancer, min_num = 15) {
    alleles_frequency_per_cancer_df %>%
        filter(cancer == !!cancer & n >= min_num) %>%
        pull(ras_allele) %>%
        unique()
}



#### ---- Mutational signature levels for all samples ---- ####

data_to_dataframe_for_clustering <- function(dat) {
    df <- dat %>%
        ungroup() %>%
        select(tumor_sample_barcode, description, contribution) %>%
        # mutate(description = paste0("sig", description)) %>%
        pivot_wider(names_from = description,
                    values_from = contribution) %>%
        as.data.frame() %>%
        column_to_rownames("tumor_sample_barcode")
    return(df)
}


dist_hclust <- function(df) {
    hc <- df %>%
        dist(method = "maximum") %>%
        hclust(method = "ward.D")
    return(hc)
}

scale_dist_hclust <- function(df) {
    scaled_df <- t(apply(df, 1, scale))
    scaled_df <- as.data.frame(scaled_df)
    colnames(scaled_df) <- colnames(df)
    scaled_df[is.na(scaled_df)] <- 0
    return(dist_hclust(scaled_df))
}

order_samples_by_hclust <- function(dat) {
    hc <- data_to_dataframe_for_clustering(dat)%>%
        scale_dist_hclust()

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

order_signatures_by_hclust <- function(dat) {
    hc <- data_to_dataframe_for_clustering(dat) %>%
        t() %>%
        scale_dist_hclust()

    sample_order_df <- tibble(
        description = hc$labels,
        signature_idx = hc$order
    )

    mod_dat <- full_join(dat, sample_order_df, by = "description")
    if(any(is.na(mod_dat))) {
        stop("Joining in `order_signatures_by_hclust()` went wrong.")
    }
    return(mod_dat)
}


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


save_barplot_distribution_per_sample_plots <- function(cancer, gg_obj, ...) {
    ggsave_wrapper(
        gg_obj,
        plot_path(GRAPHS_DIR, glue("signature-level-per-sample_{cancer}.jpeg")),
        "wide"
    )
}


mutsig_per_sample_plots <- mutsig_noartifact_df %>%
    filter(!all_zeros) %>%
    group_by(cancer, description) %>%
    filter(!all(contribution == 0)) %>%
    group_by(tumor_sample_barcode, cancer, description) %>%
    summarise(contribution = sum(contribution)) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        data = purrr::map(data, order_signatures_by_hclust),
        data = purrr::map(data, order_samples_by_hclust),
        data = purrr::map(data,
                          prep_barplot_distribution_per_sample,
                          order_signatures = TRUE),
        gg_obj = purrr::map2(data, cancer, barplot_distribution_per_sample)
    ) %>%
    pwalk(save_barplot_distribution_per_sample_plots)



mutsig_per_sample_plot <- mutsig_noartifact_df %>%
    filter(!all_zeros) %>%
    group_by(tumor_sample_barcode, cancer, description) %>%
    summarise(contribution = sum(contribution)) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        data = purrr::map(data, order_signatures_by_hclust),
        data = purrr::map(data, order_samples_by_hclust),
        data = purrr::map(data, prep_barplot_distribution_per_sample)
    ) %>%
    unnest(data) %>%
    barplot_distribution_per_sample(title = NULL) +
    facet_wrap(~cancer, scales = "free", ncol = 1)

ggsave_wrapper(
    mutsig_per_sample_plot,
    plot_path(GRAPHS_DIR, "signature-level-per-sample.jpeg"),
    "large"
)

ggsave_wrapper(
    mutsig_per_sample_plot,
    plot_path(GRAPHS_DIR, "signature-level-per-sample.svg"),
    "large"
)

saveRDS(
    mutsig_per_sample_plot,
    get_fig_proto_path("signature-level-per-sample", 2, supp = TRUE)
)


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


#### ---- Mutational signatures per KRAS allele ---- ####


# Make a barplot of the mutational signature distribution per allele.
barplot_by_allele <- function(cancer, data, ...) {
    p <- data %>%
        filter(ras_allele %in% alleles_to_plot(cancer = !!cancer)) %>%
        group_by(ras_allele, description) %>%
        summarise(contribution = sum(contribution)) %>%
        ungroup() %>%
        ggplot(aes(x = ras_allele, y = contribution)) +
        geom_col(aes(fill = description), position = "fill") +
        scale_fill_manual(
            values = mutsig_descrpt_pal,
            guide = guide_legend(
                nrow = 1,
                title.position = "left",
                label.position = "top",
                label.hjust = 0.5,
                label.vjust = -5
            )
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw(
            base_size = 9,
            base_family = "Arial"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "bottom",
            legend.key.size = unit(4, "mm"),
            legend.spacing.y = unit(0, "mm")
        ) +
        ggtitle(cancer)
    return(p)
}

# Combine the plots into a grid using 'patchwork'.
distribution_plots <- mutational_signatures_df %>%
    filter(description != "artifact" & !is_hypermutant) %>%
    mutate(
        signature = factor(signature, levels = names(mutsig_pal)),
        description = factor(description, levels = names(mutsig_descrpt_pal)),
        ras_allele = str_remove_all(ras_allele, "KRAS_"),
        ras_allele = factor(ras_allele, levels = names(short_allele_pal))
    ) %>%
    group_by(cancer) %>%
    nest() %>%
    pmap(barplot_by_allele)

combined_plots <-  wrap_plots(distribution_plots) / guide_area() +
    plot_layout(guides = 'collect', heights = c(14, 1))
ggsave_wrapper(
    combined_plots,
    plot_path(GRAPHS_DIR, "mutsig-dist_combined.svg"),
    'wide'
)



#### ---- Clock vs. non-clock ---- ####

clock_sigs <- c("1", "5")
clock_pal <- c(
    "clock" = "grey25",
    "non-clock" = "grey80"
)

plot_clock_vs_nonclock <- function(cancer, data) {
    plot_data <- data %>%
        mutate(clock_signature = ifelse(
            signature %in% !!clock_sigs, "clock", "non-clock"
        )) %>%
        group_by(tumor_sample_barcode, clock_signature) %>%
        summarise(contribution = sum(contribution)) %>%
        ungroup()


    p <- plot_data %>%
        ggplot(aes(x = clock_signature, y = contribution)) +
        geom_boxplot(
            aes(fill = clock_signature),
            alpha = 0.7,
            outlier.shape = NA,
            notch = FALSE
        ) +
        scale_fill_manual(
            values = clock_pal
        ) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
        theme_bw(
            base_size = 9,
            base_family = "Arial"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "right",
            legend.title = element_blank()
        ) +
        ggtitle(cancer)
    return(p)
}

clock_plots <- mutsig_noartifact_df %>%
    filter(!is_hypermutant) %>%
    group_by(cancer) %>%
    nest() %>%
    pmap(plot_clock_vs_nonclock)
clock_combined_plots <- wrap_plots(clock_plots) +
    plot_layout(guides = 'collect')
ggsave_wrapper(
    clock_combined_plots,
    plot_path(GRAPHS_DIR, "clock-signatures.svg"),
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
        allele = factor(allele, levels = names(short_allele_pal)),
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
