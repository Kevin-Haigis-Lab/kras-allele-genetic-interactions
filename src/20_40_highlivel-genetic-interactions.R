
#### ---- UpSetR plot for similarity of alleles ---- ####

info(logger, "Making UpSet plots for each cancer and genetic interactions.")

library(UpSetR)

# Make UpSetR plots for data nested by `cancer` and `genetic_interaction`
allele_gene_interaction_upsetr <- function(cancer, genetic_interaction, data) {

    gl <- data %>%
        group_by(allele) %>%
        summarise(genes = list(hugo_symbol)) %>%
        deframe()
    save_path <- plot_path("20_40_highlivel-genetic-interactions",
                           glue("{cancer}_{genetic_interaction}_upset.svg"))

    sizes <- c(8, 6)
    nsets <- n_distinct(names(gl))

    p <- upset(fromList(gl), nsets = nsets, nintersects = NA, order.by = "freq")
    svg(file = save_path, width = sizes[[1]], height = sizes[[2]])
    print(p)
    dev.off()
}

genetic_interaction_df %>%
    select(cancer, allele, hugo_symbol, genetic_interaction) %>%
    unique() %>%
    group_by(cancer, genetic_interaction) %>%
    nest() %>%
    pwalk(allele_gene_interaction_upsetr)



#### ---- Make high-level network images ---- ####

library(ggraph)

cache("genetic_interaction_gr",
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

set.seed()
for (CANCER in unique(genetic_interaction_df$cancer)) {
    gr_plot <- genetic_interaction_gr %E>%
        filter(cancer == !!CANCER) %N>%
        filter(centrality_degree(mode = "all") > 0) %>%
        mutate(node_label = ifelse(is_kras, name, NA),
               node_label = str_remove_all(node_label, "KRAS_")) %>%
        ggraph(layout = "stress") +
        geom_edge_link(aes(color = genetic_interaction), width = 0.1) +
        scale_edge_color_manual(values = comut_mutex_pal) +
        geom_node_point(aes(color = node_label), size = 1) +
        geom_node_text(aes(label = node_label), repel = TRUE, family = "Arial") +
        scale_color_manual(values = short_allele_pal, na.value = NA) +
        theme_graph() +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(hjust = 0.5)
        ) +
        labs(
            title = glue("Allele-specific genetic interactions in {CANCER}"),
            color = "KRAS allele",
            edge_color = "genetic\ninteraction"
        )
    save_path <- plot_path("20_40_highlivel-genetic-interactions",
                           glue("genetic_interaction_network_{CANCER}.svg"))
    ggsave_wrapper(gr_plot, save_path, size = "large")
}

