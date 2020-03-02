
GRAPHS_DIR <- "90_03_mutation-burden-distribution"
reset_graph_directory(GRAPHS_DIR)


# A simple wrapper to save figure protos from this script.
save_fig_proto_wrapper <- function(p, n) {
    proto_name <- basename(n)
    proto_name <- file_sans_ext(proto_name)
    saveRDS(
        p,
        get_fig_proto_path(proto_name, 4, supp = TRUE)
    )
}


#### ---- Find hypermutant cut-off for COAD ---- ####

upper_hypermut_bound <- cancer_coding_muts_df %>%
    filter(cancer == "COAD" & is_hypermutant) %>%
    count(tumor_sample_barcode, dataset) %>%
    group_by(dataset) %>%
    summarise(min_hypermut_muts = min(n))

lower_hypermut_bound <- cancer_coding_muts_df %>%
    filter(cancer == "COAD" & !is_hypermutant) %>%
    count(tumor_sample_barcode, dataset) %>%
    group_by(dataset) %>%
    summarise(min_hypermut_muts = max(n))



#### ---- Mut. burden distribution ---- ####
# Plots showing the distribution of the mutations per sample.

# Plot the distribution of mutational burden per sample.
# (Also saves to SVG and proto for figures.)
plot_distribution_of_mutation_count <- function(cancer,
                                                data,
                                                save_name,
                                                hypermutant_rug = FALSE) {
    p <- data %>%
        mutate(tumor_sample_barcode = fct_reorder(tumor_sample_barcode, n),
               dataset_label = paste0(dataset, "\n", target)) %>%
        ggplot(
            aes(x = tumor_sample_barcode,
                y = log10(n))
        ) +
        facet_grid(~ dataset_label, scales = "free_x") +
        geom_point(
            aes(color = log10(n)),
            size = 0.1
        ) +
        scale_color_viridis_c(
            guide = FALSE
        ) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.title.y = element_markdown(),
            panel.grid.major.x = element_blank(),
            strip.background = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        labs(
            title = as.character(
                glue("Distribution of mutations in {cancer} samples")
            ),
            x = "tumor samples",
            y = "num. mutations (*log*<sub>10</sub>-transformed)"
        )

    if (hypermutant_rug) {
        p <- p +
            geom_rug(
                aes(alpha = is_hypermutant),
                sides = "b",
                color = "black"
            ) +
            scale_alpha_manual(
                values = c("TRUE" = 0.2, "FALSE" = 0.0),
                guide = FALSE
            )
    }
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
    mutate(
        save_name = paste0(cancer, "_coding_muts_distribution.svg"),
        hypermutant_rug = cancer == "COAD"
    ) %>%
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
        scale_x_discrete(expand = c(0, 0)) +
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
            y = "fraction of mutations",
            fill = "mutation type"
        )

    save_path <- plot_path("90_03_mutation-burden-distribution", save_name)
    ggsave_wrapper(p, save_path, width = 6, height = 4)

    # Save the ggplot objects for figures.
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
