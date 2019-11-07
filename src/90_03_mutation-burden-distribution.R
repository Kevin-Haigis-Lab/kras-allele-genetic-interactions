
#### ---- Mut. burden distribution ---- ####
# Plots showing the distribution of the mutations per sample.

plot_distribution_of_mutation_count <- function(cancer, data, save_name) {
    p <- data %>%
        mutate(tumor_sample_barcode = fct_reorder(tumor_sample_barcode, n),
               dataset_label = paste0(dataset, "\n", target)) %>%
        ggplot(aes(x = tumor_sample_barcode, y = log10(n))) +
        facet_grid(~ dataset_label, scales = "free_x") +
        geom_point(aes(color = log10(n), shape = is_hypermutant), size = 0.1) +
        scale_color_viridis_c() +
        scale_shape_manual(values = c(20, 4)) +
        theme_bw(base_size = 8, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.background = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        labs(
            title = glue("Distribution of mutations in {cancer} samples"),
            x = "tumor samples",
            y = "log10( num. mutations )"
        )
    save_path <- plot_path("90_03_mutation-burden-distribution", save_name)
    ggsave_wrapper(p, save_path, "wide")
}

# Distribution of all mutations per sample.
cancer_muts_df %>%
    group_by(cancer, tumor_sample_barcode, dataset, target,
             is_hypermutant, ras_allele) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_all_muts_distribution.svg")) %>%
    pwalk(plot_distribution_of_mutation_count)

# Distribution of coding mutations per sample.
cancer_coding_muts_df %>%
    group_by(cancer, tumor_sample_barcode, dataset, target,
             is_hypermutant, ras_allele) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_coding_muts_distribution.svg")) %>%
    pwalk(plot_distribution_of_mutation_count)



#### ---- Distribution of mutation types per data source ---- ####

plot_distribution_of_mutation_type <- function(cancer, data, save_name) {
    mutation_type_hr_order <- data %>%
        group_by(mutation_type_hr) %>%
        summarise(tot = sum(n)) %>%
        ungroup() %>%
        arrange(tot) %>%
        pull(mutation_type_hr) %>%
        unlist()

    p <- data %>%
        mutate(dataset_label = paste0(dataset, "\n", target),
               mutation_type_hr = factor(mutation_type_hr, levels = mutation_type_hr_order)) %>%
        ggplot(aes(x = dataset_label, y = n)) +
        geom_col(aes(fill = mutation_type_hr), position = "fill", color = "black") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = mutation_pal) +
        theme_bw(base_size = 8, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 40, hjust = 1.0),
            panel.grid.major.x = element_blank(),
        ) +
        labs(
            title = glue("Distribution of mutations types in {cancer} samples"),
            x = "data source",
            y = "number of mutations detected",
            fill = "mutation type"
        )
    save_path <- plot_path("90_03_mutation-burden-distribution", save_name)
    ggsave_wrapper(p, save_path, "medium")
}


# Distribution of coding mutations per sample.
cancer_coding_muts_df %>%
    group_by(cancer, dataset, target, mutation_type_hr) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_mutation_types.svg")) %>%
    pwalk(plot_distribution_of_mutation_type)
