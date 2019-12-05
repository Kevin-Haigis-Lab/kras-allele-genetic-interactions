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

barplot_by_allele <- function(cancer, data, ...) {
    p <- data %>%
        filter(ras_allele %in% alleles_to_plot(cancer = !!cancer)) %>%
        group_by(ras_allele, signature) %>%
        summarise(contribution = sum(contribution)) %>%
        ungroup() %>%
        ggplot(aes(x = ras_allele, y = contribution)) +
        geom_col(aes(fill = signature), position = "fill") +
        scale_fill_manual(
            values = mutsig_pal,
            guide = guide() # TODO: 
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw(
            base_size = 8,
            base_family = "Arial"
        ) +
        theme(
            axis.title.x = element_blank()
        )
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, glue("mutsig-dist_{cancer}.svg")),
        "wide"
    )
    return(p)
}


distribution_plots <- mutational_signatures_df %>%
    mutate(
        signature = factor(signature, levels = names(mutsig_pal)),
        ras_allele = str_remove_all(ras_allele, "KRAS_")
    ) %>%
    group_by(cancer) %>%
    nest() %>%
    pmap(barplot_by_allele)


combined_plots <- wrap_plots(distribution_plots) +
    plot_layout(guides = 'collect')

ggsave_wrapper(
    combined_plots,
    plot_path(GRAPHS_DIR, "mutsig-dist_combined.svg"),
    'wide'
)
