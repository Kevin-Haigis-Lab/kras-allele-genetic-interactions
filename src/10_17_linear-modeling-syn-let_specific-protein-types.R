
# Plots of genes from specific types of genes/proteins.
# groups: transcription factors, Hubs of PPIN, cell-cycle regulators, kinases

GRAPHS_DIR <- "10_17_linear-modeling-syn-let_specific-protein-types"
reset_graph_directory(GRAPHS_DIR)

#### ---- Compile gene sets ---- ####

genesets <- list(
    transcription_factors = {
        chea_geneset_df %>%
            count(gene_set) %>%
            filter(n >= 15) %>%
            u_pull(gene_set)
    },
    kinases = {
        kea_geneset_df %>%
            count(gene_set) %>%
            filter(n >= 15) %>%
            pull(gene_set)
    },
    hubs = {
        ppiHub_geneset_df %>%
            count(gene_set) %>%
            filter(n >= 15) %>%
            pull(gene_set)
    },
    cellcycle = {
        msigdb_c2_df %>%
            filter(gene_set == "REACTOME_REGULATION_OF_MITOTIC_CELL_CYCLE") %>%
            u_pull(gene)
    }
) %>%
    enframe(name = "gene_type", value = "genes")


#### ---- Plotting protein groups ---- ####

# Make a dot-plot for all the genes in `data`.
multigene_boxplot <- function(data) {
    g <- data %>%
        group_by(hugo_symbol) %>%
        mutate(gene_effect_norm = scales::rescale(gene_effect,
                                                  to = c(-1, 1))) %>%
        ungroup() %>%
        ggplot(aes(x = hugo_symbol, y = gene_effect_norm)) +
        geom_boxplot(
            aes(fill = kras_allele, color = kras_allele),
            position = position_dodge2(),
            alpha = 0.5,
            outlier.shape = NA
        ) +
        scale_fill_manual(values = short_allele_pal) +
        scale_color_manual(values = short_allele_pal) +
        theme_bw(base_family = "Arial") +
        theme(
            axis.title.x = element_blank()
        )
    return(g)
}


# Get just the top N most variable genes.
top_n_variable_hugos <- function(df, N = 10) {
    df %>%
        group_by(hugo_symbol) %>%
        summarise(variation = var(gene_effect)) %>%
        ungroup() %>%
        top_n(10, variation) %>%
        pull(hugo_symbol)
}


# Make plots for genes of specific types:
#   - transcription factors
#   - PPIN Hubs
#   - kinases
#   - cell-cycle regulators
gene_types_plot <- function(cancer, gene_cls, hugo_symbols) {
    hugo_symbols <- unique(unlist(hugo_symbols))

    data <- depmap_modelling_df %>%
        filter_depmap_by_allele_count() %>%
        group_by(cancer) %>%
        filter(n_distinct(kras_allele) >= 3) %>%
        ungroup() %>%
        filter(cancer == !!cancer)

    purrr::pwalk(genesets, function(gene_type, genes) {
        hits <- intersect(genes, hugo_symbols)
        if (length(hits) > 0) {
            df <- data %>% filter(hugo_symbol %in% !!hits)

            if (length(hits) > 10) {
                most_var_genes <- top_n_variable_hugos(df, 10)
                df <- df %>% filter(hugo_symbol %in% !!most_var_genes)
            }

            p <- multigene_boxplot(df) +
                labs(
                    y = "effect of gene KO\n(normalized within gene)",
                    title = glue("{gene_type} in {cancer} cluster {gene_cls}")
                )
            save_path <- plot_path(
                GRAPHS_DIR,
                glue("{cancer}_{gene_type}_cluster{gene_cls}.svg")
            )
            ggsave_wrapper(p, save_path, width = 8, height = 3)
        }
    })

}


depmap_gene_clusters %>%
    group_by(cancer, gene_cls) %>%
    summarise(hugo_symbols = list(hugo_symbol)) %>%
    purrr::pwalk(gene_types_plot)



#### ---- Essential or non-essential? ---- ####
# Are any of the genes found to have varying dependency scores across KRAS
# allele listed as nonessential or essential by DepMap?

essential_genes <- essentiality_tib %>%
    mutate(label = ifelse(label != "nonessential",
                          "essential", "nonessential")) %>%
    unique() %>%
    group_by(hugo_symbol) %>%
    filter(n() == 1) %>%
    ungroup()

tbl_save_path <- table_path(
    GRAPHS_DIR,
    "essential-genes-with-allele-specific-dep.tsv"
)
inner_join(depmap_gene_clusters, essential_genes, by = "hugo_symbol") %>%
    arrange(cancer, hugo_symbol, label) %>%
    write_tsv(tbl_save_path)