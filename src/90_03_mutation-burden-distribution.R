
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



#### ---- VAF distribution ---- ####

# A function to make the labels of a ggplot a bit prettier.
pretty_rounding <- function(x) {
    case_when(x == 0 ~ "0",
              x == 1 ~ "1",
              TRUE ~ as.character(round(x, 2)))
}

# Density plot of the distribution of VAFs.
vaf_distribution_density <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    ggplot(aes(x = VAF)) +
    facet_wrap(~ cancer, scales = "free_y") +
    geom_density(aes(color = cancer)) +
    scale_color_manual(values = cancer_palette, guide = FALSE) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(
        limits = c(0, 1),
        expand = c(0, 0),
        labels = pretty_rounding
    ) +
    theme_bw(base_size = 7, base_family = "arial") +
    theme(
        strip.background = element_blank()
    )
ggsave_wrapper(
    vaf_distribution_density,
    plot_path(GRAPHS_DIR, "vaf_distribution_density.svg"),
    "small"
)



plot_mutfreq_by_avgvaf_scatter <- function(df) {
    df %>%
        ggplot(aes(x = avg_vaf, y = mut_freq)) +
        facet_wrap(~ cancer, scales = "free") +
        geom_point(size = 0.3, alpha = 0.7, color = "grey25") +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 1.2, family = "arial", color = "blue",
            segment.color = "grey40", segment.size = 0.2, segment.alpha = 0.5,
            seed = 0,
        ) +
        scale_x_continuous(
            limits = c(0, 1),
            expand = c(0, 0),
            labels = pretty_rounding
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.02)),
            labels = pretty_rounding
        ) +
        theme_bw(base_size = 7, base_family = "arial") +
        theme(
            strip.background = element_blank()
        ) +
        labs(
            x = "average VAF",
            y = "mutation frequency"
        )
}

# Average VAF of mutations in a gene vs. mutational frequency of the gene.
vaf_mutfreq_scatter <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    group_by(cancer) %>%
    mutate(num_cancer_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer, hugo_symbol) %>%
    summarise(
        avg_vaf = mean(VAF),
        mut_freq = n_distinct(tumor_sample_barcode) / unique(num_cancer_samples)
    ) %>%
    group_by(cancer) %>%
    arrange(-mut_freq, -avg_vaf) %>%
    mutate(label = ifelse(1:n() < 20, hugo_symbol, NA)) %>%
    ungroup() %>%
    plot_mutfreq_by_avgvaf_scatter()

ggsave_wrapper(
    vaf_mutfreq_scatter,
    plot_path(GRAPHS_DIR, "vaf_mutfreq_scatter.jpeg"),
    "medium"
)


# A simple data frame of the genes with genetic interactions per cancer type.
genetic_interaction_genes <- genetic_interaction_df %>%
    select(cancer, hugo_symbol) %>%
    unique() %>%
    bind_rows(
        tibble(cancer = names(cancer_palette), hugo_symbol = "KRAS")
    )


# Average VAF of mutations in a gene vs. mutational frequency of the gene, but
# only for genes with a comutation interaction
vaf_mutfreq_comuts_scatter <- cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    inner_join(genetic_interaction_genes, by = c("cancer", "hugo_symbol")) %>%
    group_by(cancer) %>%
    mutate(num_cancer_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer, hugo_symbol) %>%
    summarise(
        avg_vaf = mean(VAF),
        mut_freq = n_distinct(tumor_sample_barcode) / unique(num_cancer_samples)
    ) %>%
    group_by(cancer) %>%
    arrange(-mut_freq, -avg_vaf) %>%
    mutate(label = ifelse(1:n() < 20, hugo_symbol, NA)) %>%
    ungroup() %>%
    plot_mutfreq_by_avgvaf_scatter()

ggsave_wrapper(
    vaf_mutfreq_comuts_scatter,
    plot_path(GRAPHS_DIR, "vaf_mutfreq_comuts_scatter.svg"),
    "medium"
)


# Simple rounding alias.
r3 <- function(x) round(x, 3)

# Summary statistics on VAF values per cancer type.
cancer_coding_muts_df %>%
    filter(cancer != "SKCM" & !is.na(VAF)) %>%
    group_by(cancer) %>%
    summarise(
        min_vaf = r3(min(VAF)),
        q25_vaf = r3(quantile(VAF, 0.25)),
        avg_vaf = r3(mean(VAF)),
        mid_vaf = r3(median(VAF)),
        q75_vaf = r3(quantile(VAF, 0.75)),
        max_vaf = r3(max(VAF)),
        stddev_vaf = r3(sd(VAF))
    ) %>%
    ungroup() %>%
    knitr::kable()
