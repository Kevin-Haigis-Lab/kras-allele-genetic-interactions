



all_transcription_factors <- unique(c(
    tf2dna_tfs,
    unique(chea_geneset_df$gene_set)
))


# Make a dot-plot for all the genes in `data`.
multigene_dotplot <- function(data) {
    g <- data %>%
        group_by(hugo_symbol) %>%
        mutate(gene_effect_norm = scales::rescale(gene_effect, to = c(-1, 1))) %>%
        ungroup() %>%
        ggplot(aes(x = hugo_symbol, y = gene_effect_norm)) +
        geom_boxplot(aes(fill = allele), outlier.shape = NA) +
        theme_bw()
    return(g)
}



# Make plots for genes of specific types:
#   - transcription factors
#   - PPIN Hubs
#   - kinases
#   - cell-cycle regulators
gene_types_plot <- function(cancer, gene_cls, hugo_symbols) {
    hugo_symbols <- unique(unlist(hugo_symbols))

    data <- model_data %>%
        filter(cancer == !!cancer)

    if (any(hugo_symbols %in% all_transcription_factors)) {
        data %>%
            filter(
                hugo_symbol %in% !!hugo_symbols &
                hugo_symbol %in% all_transcription_factors
            ) %>%
            multigene_dotplot()
    }





}



depmap_gene_clusters %>%
    group_by(cancer, gene_cls) %>%
    summarise(hugo_symbols = list(hugo_symbol)) %>%
    purrr::pwalk(gene_types_plot)