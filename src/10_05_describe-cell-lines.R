# Some descriptive plots for the cell lines used in DepMap project

library(gtable)
library(gridExtra)


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

# A place-holder for plotting.
empty_grob <- ggplotGrob(
    model_data  %>%
    sample_n(10) %>%
    ggplot(aes(x = dep_map_id, y = hugo_symbol)) +
    geom_tile(fill = NA, color = NA) +
    theme_void()
)


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
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
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
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
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

    ## The rest are concerned with plotting the above pieces in one plot.

    col_bar_grob <- ggplotGrob(col_bar_plot)
    mut_tile_grob <- ggplotGrob(mut_tile_plot)
    row_bar_grob <- ggplotGrob(row_bar_plot)
    kras_allele_grob <- ggplotGrob(kras_allele_plot)

    layout_matrix <- rbind(
        c(1, NA),
        c(2, 3),
        c(4, NA)
    )

    col_one_width <- grid::unit.pmax(col_bar_grob$widths[2:5], mut_tile_grob$widths[2:5], kras_allele_grob$widths[2:5])
    col_bar_grob$widths[2:5] <- as.list(col_one_width)
    mut_tile_grob$widths[2:5] <- as.list(col_one_width)
    kras_allele_grob$widths[2:5] <- as.list(col_one_width)

    row_two_height <- grid::unit.pmin(mut_tile_grob$heights, row_bar_grob$heights)
    mut_tile_grob$heights <- as.list(row_two_height)
    row_bar_grob$heights <- as.list(row_two_height)

    mut_tile_grid_plot <- arrangeGrob(
        grobs = list(col_bar_grob, mut_tile_grob, row_bar_grob, kras_allele_grob),
        layout_matrix = layout_matrix,
        heights = c(1, 5, ifelse(CANCER == "PAAD", 0.8, 0.6)),
        widths = c(8, 1)
    )


    ## Saving the plot
    
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

    ggsave_wrapper(mut_tile_grid_plot, save_path,
                   width = save_width, height = save_height)
}
