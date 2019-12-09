
# Analyze the comutation interaction networks at a high level.

GRAPHS_DIR <- "20_40_highlivel-genetic-interactions"
reset_graph_directory(GRAPHS_DIR)

#### ---- UpSetR plot for similarity of alleles ---- ####

# Make UpSetR plots for data nested by `cancer` and `genetic_interaction`
allele_gene_interaction_upsetr <- function(cancer, genetic_interaction, data) {

    gl <- data %>%
        group_by(allele) %>%
        summarise(genes = list(hugo_symbol)) %>%
        deframe()
    save_path <- plot_path(GRAPHS_DIR,
                           glue("{cancer}_{genetic_interaction}_upset.svg"))

    sizes <- c(8, 6)
    nsets <- n_distinct(names(gl))

    p <- upset(fromList(gl), nsets = nsets, nintersects = NA, order.by = "freq")
    svg(file = save_path, width = sizes[[1]], height = sizes[[2]])
    print(p)
    dev.off()
}

caps <- as.list(capabilities())
if (caps$X11) {
    info(logger, "Making UpSet plots for each cancer and genetic interactions.")

    library(UpSetR)

    genetic_interaction_df %>%
        select(cancer, allele, hugo_symbol, genetic_interaction) %>%
        unique() %>%
        group_by(cancer, genetic_interaction) %>%
        nest() %>%
        pwalk(allele_gene_interaction_upsetr)
}





#### ---- Make high-level network images ---- ####

library(ggraph)

ProjectTemplate::cache("genetic_interaction_gr",
                       depends = "genetic_interaction_df",
{
    genetic_interaction_gr <- genetic_interaction_df %>%
        select(hugo_symbol, kras_allele, cancer, p_val, genetic_interaction) %>%
        unique() %>%
        mutate(from = kras_allele, to = hugo_symbol) %>%
        filter(str_replace_us(kras_allele) %in% names(allele_palette)) %>%
        as_tbl_graph() %N>%
        mutate(is_kras = str_detect(name, "KRAS_"))
    return(genetic_interaction_gr)
})


goi_genes <- unique(genes_of_interest_df$hugo_symbol)

# Plots the genetic interaction networks.
make_network_plot <- function(gr) {
    p <- gr %>%
        ggraph(layout = "nicely") +
        geom_edge_link(
            aes(color = genetic_interaction),
            width = 0.3
        ) +
        scale_edge_color_manual(values = comut_mutex_pal) +
        geom_node_point(
            aes(color = node_color,
                size = node_size),
        ) +
        geom_node_text(
            aes(label = node_label,
                size = label_size),
            repel = TRUE,
            family = "Arial"
        ) +
        scale_color_manual(
            values = short_allele_pal,
            na.value = NA
        ) +
        scale_size_identity() +
        theme_graph(
            base_family = "Arial",
            base_size = 8
        ) +
        theme(
            plot.title = element_blank(),
            legend.position = "left"
        )
    return(p)
}


for (CANCER in sort(unique(genetic_interaction_df$cancer))) {
    set.seed(0)
    gr_plot <- genetic_interaction_gr %E>%
        filter(cancer == !!CANCER) %N>%
        filter(centrality_degree(mode = "all") > 0) %>%
        mutate(node_label = ifelse(is_kras, name, NA),
               node_label = str_remove_all(node_label, "KRAS_"),
               node_color = node_label,
               node_size = 2,
               label_size = 5) %>%
        make_network_plot() +
        labs(
            edge_color = "interaction"
        )
    save_path <- plot_path(GRAPHS_DIR,
                           glue("genetic_interaction_network_{CANCER}.svg"))
    ggsave_wrapper(gr_plot, save_path, size = "medium")

    saveRDS(
        gr_plot,
        get_fig_proto_path(glue("genetic_interaction_network_{CANCER}.svg"), 2)
    )
}



#### ---- Correlate number of interactors KRAS allele frequency ---- ####

interaction_counts_data <- genetic_interaction_df %>%
    mutate(allele_freq = ifelse(
        genetic_interaction == "comutation",
        allele_freq,
        num_samples_per_cancer_allele / num_samples_per_cancer)) %>%
    group_by(cancer, allele, allele_freq, genetic_interaction) %>%
    count() %>%
    ungroup()

interaction_counts_plot <- interaction_counts_data %>%
    filter(n > 0) %>%
    ggplot(aes(x = allele_freq, y = n)) +
    facet_wrap(~ cancer, scales = "free") +
    geom_point(aes(color = allele, shape = genetic_interaction)) +
    scale_color_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(0, 2)) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(breaks = integer_breaks(),
                       limits = c(0, NA)) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank()
    ) +
    labs(
        title = "Association of allele frequency and number of interactors",
        x = "KRAS allele frequency",
        y = "number of genetic interactors",
        color = "allele",
        shape = "interaction"
    )
save_path <- plot_path(GRAPHS_DIR,
                       "corr_allele-freq_num-interactors.svg")
ggsave_wrapper(interaction_counts_plot, save_path, "medium")
