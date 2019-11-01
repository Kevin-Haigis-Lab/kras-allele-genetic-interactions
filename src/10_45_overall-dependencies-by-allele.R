# The overall distribution of dependencies of cell lines by KRAS allele.

library(ggridges)

#### ---- Ridges and Quantiles ---- ####


calculate_quantiles <- function(x, n = 5) {
    quantile(x, probs = seq(0, 1, 1/n)) %>%
        enframe(name = "percentile", value = "value") %>%
        mutate(percentile = factor(percentile, levels = percentile)) %>%
        list()
}



for (CANCER in unique(model_data$cancer)) {
    if (CANCER == "MM") next
    p1 <- model_data %>%
        filter(cancer == !!CANCER & !is.na(gene_effect)) %>%
        group_by(allele, hugo_symbol) %>%
        summarise(avg_gene_effect = mean(gene_effect)) %>%
        ungroup() %>%
        filter(avg_gene_effect < 0.5) %>%
        ggplot(aes(x = avg_gene_effect, y = allele)) +
        geom_density_ridges(aes(fill = allele), color = "white", alpha = 0.7) +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        scale_y_discrete(expand = expand_scale(mult = c(0.02, 0.02))) +
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

    p2 <- model_data %>%
        filter(cancer == !!CANCER & !is.na(gene_effect)) %>%
        group_by(allele) %>%
        summarise(data = calculate_quantiles(gene_effect, 10)) %>%
        ungroup() %>%
        unnest(data) %>%
        ggplot(aes(x = percentile, y = value)) +
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        geom_line(aes(color = allele, group = allele), alpha = 0.7) +
        geom_point(aes(color = allele), alpha = 0.7) +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        theme_bw(base_family = "Arial", base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = c(0.1, 0.8),
            legend.box.background = element_rect(color = "grey25", fill = "white"),
            legend.box.margin = margin(0, 0.2, 0, 0, "lines"),
            legend.title = element_blank()
        ) +
        labs(
            x = "percentile",
            y = "depletion effect value",
            title = glue("Quantiles of DepMap scores per KRAS allele in {CANCER}")
        )

    cp <- cowplot::plot_grid(p1, p2, nrow = 1, align = "h")

    save_path <- plot_path("10_45_overall-dependencies-by-allele",
                           glue("dependency-overview_{CANCER}.svg"))
    cowplot::save_plot(save_path, cp, base_width = 10, base_height = 4)
}
