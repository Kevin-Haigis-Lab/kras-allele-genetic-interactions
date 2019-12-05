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
