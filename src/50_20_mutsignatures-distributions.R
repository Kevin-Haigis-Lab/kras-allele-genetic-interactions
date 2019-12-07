# Plot the distribution of mutational signatures by allele.

library(patchwork)

GRAPHS_DIR <- "50_20_mutsignatures-distributions"
reset_graph_directory(GRAPHS_DIR)

alleles_frequency_per_cancer_df <- mutational_signatures_df %>%
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
    filter(description != "artifact") %>%
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
