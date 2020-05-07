# The overall distribution of dependencies of cell lines by KRAS allele.

GRAPHS_DIR <- "10_45_overall-dependencies-by-allele"
reset_graph_directory(GRAPHS_DIR)

library(ggridges)

#### ---- Ridges and Quantiles ---- ####


calculate_quantiles <- function(x, n = 5) {
    quantile(x, probs = seq(0, 1, 1/n)) %>%
        enframe(name = "percentile", value = "value") %>%
        mutate(percentile = factor(percentile, levels = percentile)) %>%
        list()
}



for (CANCER in unique(depmap_modelling_df$cancer)) {
    p1 <- depmap_modelling_df %>%
        filter_depmap_by_allele_count() %>%
        filter(cancer == !!CANCER & !is.na(gene_effect)) %>%
        group_by(kras_allele, hugo_symbol) %>%
        summarise(avg_gene_effect = mean(gene_effect)) %>%
        ungroup() %>%
        filter(avg_gene_effect < 0.5) %>%
        ggplot(aes(x = avg_gene_effect, y = kras_allele)) +
        geom_density_ridges(aes(fill = kras_allele),
                            color = "white", alpha = 0.7) +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
        theme_bw(base_family = "Arial", base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
        ) +
        labs(
            x = "average depletion effect",
            y = "KRAS allele",
            title = glue("Distribution of dependencies in {CANCER} cell lines")
        )

    p2 <- depmap_modelling_df %>%
        filter_depmap_by_allele_count() %>%
        filter(cancer == !!CANCER & !is.na(gene_effect)) %>%
        group_by(kras_allele) %>%
        summarise(data = calculate_quantiles(gene_effect, 10)) %>%
        ungroup() %>%
        unnest(data) %>%
        ggplot(aes(x = percentile, y = value)) +
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        geom_line(aes(color = kras_allele, group = kras_allele), alpha = 0.7) +
        geom_point(aes(color = kras_allele), alpha = 0.7) +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        theme_bw(base_family = "Arial", base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = c(0.1, 0.8),
            legend.box.background = element_rect(color = "grey25",
                                                 fill = "white"),
            legend.box.margin = margin(0, 0.2, 0, 0, "lines"),
            legend.title = element_blank()
        ) +
        labs(
            x = "percentile",
            y = "depletion effect value",
            title = glue("Quantiles of DepMap scores per KRAS allele in {CANCER}")
        )

    patch <- p1 + p2

    ggsave_wrapper(
        patch,
        plot_path(GRAPHS_DIR, glue("dependency-overview_{CANCER}.svg")),
        width = 10, height = 4
    )
}
