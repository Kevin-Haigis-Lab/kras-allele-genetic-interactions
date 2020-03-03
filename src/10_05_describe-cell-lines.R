# Some descriptive plots for the cell lines used in DepMap project


GRAPHS_DIR <- "10_05_describe-cell-lines"
reset_graph_directory(GRAPHS_DIR)


#### ---- Shared mutations ---- ####
# Mutations shared by cell lines.

# Genes to include in the plot.
cgc_genes <-  cosmic_cgc_df %>%
    filter(tier == 1) %>%
    u_pull(hugo_symbol)


# The cancer and KRAS allele for each cell line.
ccle_cancers <- model_data %>%
    select(dep_map_id, cancer, allele) %>%
    unique()


# Mutation data to plot.
mut_df <- ccle_mutations_coding %>%
    filter(dep_map_id %in% !!model_data$dep_map_id) %>%
    filter(hugo_symbol %in% !!cgc_genes) %>%
    inner_join(ccle_cancers, by = "dep_map_id") %>%
    filter(cancer != "MM" & hugo_symbol != "KRAS")


# Create a plot for each cancer.
for (CANCER in unique(mut_df$cancer)) {

    # Text color for a dot representing CGC hotspot mutations.
    label_color <- ifelse(CANCER == "COAD", "black", "white")

    # Data to plot for the cancer.
    # `x` and `y` are factors to control the order of plotting.
    mut_cancer_data <- mut_df %>%
        filter(cancer == !!CANCER) %>%
        complete(dep_map_id, hugo_symbol) %>%
        arrange(allele) %>%
        mutate(x = fct_inorder(dep_map_id),
               y = factor(hugo_symbol, rev(sort(unique(hugo_symbol)))),
               mut_label = ifelse(is_cosmic_hotspot, "•", ""))

    # Bar plot to show the margin for the columns.
    col_bar_plot <- mut_cancer_data %>%
        filter(!is.na(cancer)) %>%
        count(x) %>%
        ggplot(aes(x = x, y = n)) +
        geom_col(aes(alpha = n), fill = "dodgerblue4") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
        theme_classic(base_family = "Arial", base_size = 5) +
        theme(
            legend.position = 'none',
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        labs(
            title = glue("CGC mutations in {CANCER} cell lines")
        )

    # Tile plot for which cell lines have a mutation in each CGC gene.
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
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        labs(
            y = "CGC genes (dot = hotspot mutation)"
        )

    # Bar plot to show the margin for the rows.
    row_bar_plot <- mut_cancer_data %>%
        filter(!is.na(cancer)) %>%
        count(y) %>%
        ggplot(aes(x = y, y = n)) +
        geom_col(aes(alpha = n), fill = "dodgerblue4") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
        theme_classic(base_family = "Arial", base_size = 5) +
        theme(
            legend.position = 'none',
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        coord_flip()

    # A strip showing KRAS allele along the top
    kras_allele_plot <- mut_cancer_data %>%
        ggplot(aes(x = x, y = "KRAS")) +
        geom_tile(aes(fill = allele), color = "black") +
        scale_fill_manual(values = short_allele_pal, guide = FALSE) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        theme_bw(base_family = "Arial", base_size = 5) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_text(angle = 45, hjust = 1.0),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        labs(
            y = "cell line (DepMap ID)"
        )

    # Saving the plot
    save_path <- plot_path(GRAPHS_DIR, glue("gcg-mutation-tile_{CANCER}.svg"))
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

    design <- "
        11111##
        2222233
        2222233
        2222233
        2222233
        2222233
        2222233
        44444##
    "

    widths <- c(10, 2)
    mut_tile_patch <- (
            ((col_bar_plot | plot_spacer()) + plot_layout(widths = widths)) /
            ((mut_tile_plot | row_bar_plot) + plot_layout(widths = widths)) /
            ((kras_allele_plot | plot_spacer()) + plot_layout(widths = widths))
        ) +
        plot_layout(heights = c(5, 20, 1))

    ggsave_wrapper(mut_tile_patch, save_path,
                   width = save_width, height = save_height)
}


#### ---- Effect of KRAS KO ---- ####

kras_gene_effect_boxplot <- model1_tib %>%
    filter(hugo_symbol == "KRAS") %>%
    select(hugo_symbol, cancer, data) %>%
    unnest(data)  %>%
    ggplot(aes(x = allele, y = gene_effect, color = allele)) +
    facet_grid(~ cancer, scales = "free_x") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.8) +
    scale_color_manual(values = short_allele_pal) +
    theme_bw(base_size = 7, base_family = "arial") +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        axis.title.y = element_markdown(),
        axis.title.x = element_blank()
    ) +
    labs(y = "depletion effect of *KRAS* KO")
ggsave_wrapper(
    kras_gene_effect_boxplot,
    plot_path(GRAPHS_DIR, "kras-gene-effect_all-cancers_boxpolot.svg"),
    "wide"
)
