# A file for making a few miscellaneous plots for presentations, etc.

GRAPHS_DIR <- "90_34_random-description-plots"
reset_graph_directory(GRAPHS_DIR)



#### ---- DepMap value distribution ---- ####

depmap_dist_plot_data <- depmap_modelling_df %>%
    filter(cancer == "PAAD") %>%
    group_by(hugo_symbol) %>%
    summarize(gene_effect = median(gene_effect, na.rm = TRUE)) %>%
    ungroup()

depmap_dist_plot_data2 <- depmap_dist_plot_data %>%
    inner_join(
        essentiality_tib %>%
            select(hugo_symbol, label) %>%
            filter(label %in% c("achilles_essential", "nonessential")) %>%
            group_by(hugo_symbol) %>%
            filter(n() == 1 | label == "nonessential") %>%
            ungroup(),
        by = "hugo_symbol"
    ) %>%
    mutate(label = case_when(
        label == "achilles_essential" ~ "essential genes",
        label == "nonessential" ~ "nonessential genes",
        TRUE ~ NA_character_
    ))


essential_pal <- c(
    "essential genes" = "#1B9E77",
    "nonessential genes" = "#D95F02"
)


depmap_dist_plot <- depmap_dist_plot_data %>%
    ggplot(aes(x = gene_effect)) +
    geom_vline(xintercept = c(0, -1), lty = 2, color = "grey50", size = 1) +
    geom_density(aes(color = label, fill = label), alpha = 0.5, size = 0,
                 data = depmap_dist_plot_data2) +
    geom_density(color = "grey20", fill = NA, size = 1.1) +
    annotate("text",
             label = names(essential_pal), color = essential_pal,
             x = c(-2, -0.5),
             y = c(0.6, 2.5),
             family = "Arial", size = 6, fontface = "bold"
         ) +
    scale_x_continuous(limits = c(NA, 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = essential_pal, guide = FALSE) +
    scale_color_manual(values = essential_pal, guide = FALSE) +
    theme_minimal(base_size  = 14, base_family = "Arial") +
    theme(
        legend.position = c(0.1, 0.2),
        legend.title = element_blank()
    ) +
    labs(x = "median dependency score per gene",
         y = "density")
ggsave_wrapper(
    depmap_dist_plot,
    plot_path(GRAPHS_DIR, "depmap-dist-plot.jpeg"),
    "wide"
)


