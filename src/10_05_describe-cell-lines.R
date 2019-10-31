# Some descriptive plots for the cell lines used in DepMap project

#### ---- Shared mutations ---- ####
# Mutations shared by cell lines.

cgc_genes <-  cosmic_cgc_df %>%
    filter(tier == 1) %>%
    u_pull(hugo_symbol)

ccle_cancers <- model_data %>%
    select(dep_map_id, cancer, allele) %>%
    unique()


mut_df <- ccle_mutations_coding %>%
    filter(dep_map_id %in% !!model_data$dep_map_id) %>%
    filter(hugo_symbol %in% !!cgc_genes) %>%
    inner_join(ccle_cancers, by = "dep_map_id") %>%
    filter(cancer != "MM" & hugo_symbol != "KRAS")

for (CANCER in unique(mut_df$cancer)) {
    label_color <- ifelse(CANCER == "COAD", "black", "white")

    # data to plot
    mut_cancer_data <- mut_df %>%
        filter(cancer == !!CANCER) %>%
        complete(dep_map_id, hugo_symbol) %>%
        arrange(allele) %>%
        mutate(x = fct_inorder(dep_map_id),
               y = factor(hugo_symbol, rev(sort(unique(hugo_symbol)))),
               mut_label = ifelse(is_cosmic_hotspot, "•", ""))

    # A strip showing KRAS allele along the top
    kras_allele_plot <- mut_cancer_data %>%
        ggplot(aes(x = x, y = "KRAS")) +
        geom_tile(aes(fill = allele), color = "black") +
        scale_fill_manual(values = short_allele_pal, guide = FALSE) +
        theme_bw(base_family = "Arial", base_size = 5) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()
        ) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        labs(title = glue("CGC mutations in {CANCER} cell lines"))


    mut_tile_plot <- mut_cancer_data %>%
        ggplot(aes(x = x, y = y)) +
        geom_tile(aes(fill = cancer), color = "black") +
        geom_text(aes(label = mut_label),
                  color = label_color,
                  size = 3) +
        scale_fill_manual(
            values = cancer_palette, na.value = "white",
            guide = FALSE) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        theme_bw(base_family = "Arial", base_size = 5) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_text(angle = 45, hjust = 1.0),
            panel.grid.major = element_blank()
        ) +
        labs(
            x = "cell line (DepMap ID)",
            y = "CGC genes (dot = hotspot mutation)"
        )

    save_path <- plot_path("10_05_describe-cell-lines",
                           glue("gcg-mutation-tile_{CANCER}.svg"))
    save_width <- case_when(
        CANCER == "LUAD" ~ 8,
        CANCER == "PAAD" ~ 4,
        CANCER == "COAD" ~ 6
    )
    save_height <- case_when(
        CANCER == "LUAD" ~ 4,
        CANCER == "PAAD" ~ 3,
        CANCER == "COAD" ~ 4
    )

    mut_tile_cowplot <- cowplot::plot_grid(
        kras_allele_plot, mut_tile_plot,
        align = "v",
        ncol = 1,
        rel_heights = c(1, ifelse(CANCER == "PAAD", 10, 13))
    )

    cowplot::save_plot(save_path, mut_tile_cowplot,
                       base_width = save_width, base_height = save_height)
}
