# Effects of KO-ing NF1 in COAD cell lines.

p <- depmap_modelling_df %>%
    filter_depmap_by_allele_count() %>%
    group_by(cancer) %>%
    filter(n_distinct(kras_allele) >= 3) %>%
    ungroup() %>%
    filter(cancer == "COAD" & hugo_symbol == "NF1") %>%
    unique() %>%
    mutate(dep_map_id = fct_reorder(dep_map_id, -gene_effect)) %>%
    ggplot() +
    geom_col(
        aes(
            x = dep_map_id,
            y = gene_effect,
            fill = kras_allele,
            alpha = is_mutated
    )) +
    scale_fill_manual(
        values = short_allele_pal,
        guide = guide_legend(title.position = "top")
    ) +
    scale_alpha_manual(
        values = c(0.5, 1),
        guide = guide_legend(title.position = "top")
    ) +
    theme_bw() +
    theme(
        text = element_text("arial"),
        axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "bottom"
    ) +
    labs(
        title = "The effects of knocking out NF1 in COAD",
        x = "cell line (DepMap ID)",
        y = "effect on growth",
        alpha = "NFI is altered",
        fill = "KRAS allele"
    )
save_path <- plot_path("90_10_NF1-depletion-COAD",
                       "NF1-depletion-COAD-barplot.svg")
ggsave_wrapper(p, save_path, size = "medium")

