# A few plots describing the results of removing unexpressed genes.

GRAPHS_DIR <- "05_10_expression-filter-plots"
TABLES_DIR <- "05_10_expression-filter-plots"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(TABLES_DIR)


unexprsd_genes_tib <- tibble::enframe(confidently_unexpressed_genes) %>%
    unnest(value) %>%
    dplyr::rename(cancer = name,
                  hugo_symbol = value) %>%
    mutate(cancer = str_to_upper(cancer)) %>%
    filter(cancer != "SKCM")


mutation_count_tib <- cancer_coding_muts_df %>%
    group_by(cancer) %>%
    mutate(num_samples = n_distinct(tumor_sample_barcode)) %>%
    group_by(cancer, hugo_symbol, num_samples) %>%
    summarise(num_mutations = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    mutate(freq_mutations = num_mutations / num_samples)


unexprsd_genes_tib %<>%
    left_join(mutation_count_tib, by = c("cancer", "hugo_symbol")) %>%
    mutate(
        num_mutations = replace_na_zero(num_mutations),
        freq_mutations = replace_na_zero(freq_mutations)
    )


mut_freq_hist <- unexprsd_genes_tib %>%
    ggplot(aes(x = freq_mutations)) +
    facet_wrap(~ cancer, scales = "free", ncol = 2) +
    geom_histogram(
        aes(fill = cancer, color = cancer),
        bins = 20,
        alpha = 0.5
    ) +
    scale_fill_manual(values = cancer_palette) +
    scale_color_manual(values = cancer_palette) +
    scale_x_continuous(
        expand = expand_scale(mult = c(0, 0.02)),
    ) +
    scale_y_continuous(
        expand = expand_scale(mult = c(0, 0.02)),
        breaks = c(10, 100, 1000, 10000, 40000)
    ) +
    coord_trans(y = my_trans_log10) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        strip.background = element_blank(),
        legend.position = "none"
    ) +
    labs(
        x = "frequency of mutation",
        y = "count"
    )

ggsave_wrapper(
    mut_freq_hist,
    plot_path(GRAPHS_DIR, "mutation-frequency-hist.svg"),
    width = 6, height = 4
)
saveRDS(
    mut_freq_hist,
    get_fig_proto_path("mutation-frequency-hist", 3, supp = TRUE)
)


#### ---- Data to TSV ---- ####

# All data to a single TSV for each source.
write_tsv(
    gtex_summary_expr,
    table_path(TABLES_DIR, "GTEx_summary.tsv")
)

write_tsv(
    hpa_expr,
    table_path(TABLES_DIR, "HPA_summary.tsv")
)

write_tsv(
    normal_tissue_expr,
    table_path(TABLES_DIR, "GTEx-and-HPA_summary.tsv")
)

write_tsv(
    cancer_rna_tib,
    table_path(TABLES_DIR, "cancer-RNA_summary.tsv")
)

# Summary of number of samples for GTEx and cancer samples.
# HPA only provides summary values per tissue.
gtex_summary_expr %>%
    group_by(tissue) %>%
    summarise(num_samples = max(GTEx_num_samples)) %>%
    filter(tissue != "skin") %>%
    ungroup() %>%
    arrange(tissue) %>%
    knitr::kable()

cancer_rna_tib %>%
    mutate(cancer = str_to_upper(cancer)) %>%
    group_by(cancer) %>%
    summarise(num_samples = max(num_samples)) %>%
    filter(cancer != "SKCM") %>%
    ungroup() %>%
    arrange(cancer) %>%
    knitr::kable()