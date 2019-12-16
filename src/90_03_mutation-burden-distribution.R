
GRAPHS_DIR <- "90_03_mutation-burden-distribution"
reset_graph_directory(GRAPHS_DIR)


# A simple wrapper to save figure protos from this script.
save_fig_proto_wrapper <- function(p, n) {
    proto_name <- basename(n)
    proto_name <- tools::file_path_sans_ext(proto_name)
    saveRDS(
        p,
        get_fig_proto_path(proto_name, 4, supp = TRUE)
    )
}


# Rename some specific data sets.
rename_datasets <- function(ds) {
    str_remove_all(ds, "coadread_|coad_|luad_|mm_|paad_") %>%
        str_remove_all("pan_can_atlas") %>%
        str_replace("__", "_")
}
rename_datasets <- memoise::memoise(rename_datasets)



#### ---- Mut. burden distribution ---- ####
# Plots showing the distribution of the mutations per sample.

plot_distribution_of_mutation_count <- function(cancer, data, save_name) {
    p <- data %>%
        mutate(tumor_sample_barcode = fct_reorder(tumor_sample_barcode, n),
               dataset_label = paste0(dataset, "\n", target)) %>%
        ggplot(
            aes(x = tumor_sample_barcode,
                y = log10(n))
        ) +
        facet_grid(~ dataset_label, scales = "free_x") +
        geom_point(
            aes(color = log10(n),
                shape = is_hypermutant),
            size = 0.1
        ) +
        scale_color_viridis_c(
            guide = FALSE
        ) +
        scale_shape_manual(
            values = c(20, 4),
            guide = FALSE
        ) +
        theme_bw(base_size = 7, base_family = "Arial") +
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
    save_path <- plot_path(GRAPHS_DIR, save_name)
    ggsave_wrapper(p, save_path, width = 8, height = 2.5)

    # Save for use in a figure.
    save_fig_proto_wrapper(p, save_path)
}

# Distribution of all mutations per sample.
cancer_muts_df %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    group_by(cancer, tumor_sample_barcode, dataset, target,
             is_hypermutant, ras_allele) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_all_muts_distribution.svg")) %>%
    pwalk(plot_distribution_of_mutation_count)

# Distribution of coding mutations per sample.
cancer_coding_muts_df %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    group_by(cancer, tumor_sample_barcode, dataset, target,
             is_hypermutant, ras_allele) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_coding_muts_distribution.svg")) %>%
    pwalk(plot_distribution_of_mutation_count)



#### ---- Distribution of mutation types per data source ---- ####

plot_distribution_of_mutation_type <- function(cancer, data, save_name) {
    # Get the ordering for the mutation types.
    mutation_type_hr_order <- data %>%
        group_by(mutation_type_hr) %>%
        summarise(tot = sum(n)) %>%
        ungroup() %>%
        arrange(tot) %>%
        pull(mutation_type_hr) %>%
        unlist()

    p <- data %>%
        mutate(
            dataset_label = paste0(dataset, "\n", target),
            mutation_type_hr = factor(mutation_type_hr,
                                      levels = mutation_type_hr_order)
        ) %>%
        ggplot(aes(x = dataset_label, y = n)) +
        geom_col(
            aes(fill = mutation_type_hr),
            position = "fill"
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(
            values = mutation_pal,
            guide = guide_legend(
                nrow = 2
            )
        ) +
        theme_bw(base_size = 8, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 40, hjust = 1.0),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        labs(
            y = "number of mutations detected",
            fill = "mutation type"
        )

    save_path <- plot_path("90_03_mutation-burden-distribution", save_name)
    ggsave_wrapper(p, save_path, width = 6, height = 4)

    save_fig_proto_wrapper(p, save_path)

    return(p)
}


# Distribution of coding mutations per sample.
mut_type_bars <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM") %>%
    mutate(dataset = rename_datasets(dataset)) %>%
    group_by(cancer, dataset, target, mutation_type_hr) %>%
    count() %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(save_name = paste0(cancer, "_mutation_types.svg")) %>%
    pwalk(plot_distribution_of_mutation_type)
