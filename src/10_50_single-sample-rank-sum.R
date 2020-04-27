# A rank-sum test within each cell line for genes in relevant gene sets.

library(gtable)
library(gridExtra)


#### ---- Prepare gene-sets ---- ####

# Compress a gene set tibble into a single row per gene set.
compress_gene_set_tibble <- function(gs_tib, gene_col) {
    gene_col <- rlang::enquo(gene_col)
    gs_tib %>%
        group_by(gene_set) %>%
        summarise(genes = list(!!gene_col))
}

# All gene sets in one tibble with one row per gene set.
gene_sets_df <- bind_rows(
    compress_gene_set_tibble(kegg_geneset_df, hugo_symbol),
    compress_gene_set_tibble(kea_geneset_df, gene),
    compress_gene_set_tibble(filter(msigdb_c2_df, str_detect(gene_set, "REACTOME")))
)


# Wilcoxon rank-sum test for the genes in a gene set in a single sample.
ranksum_enrichment <- function(df, gene_set, ...) {
    x <- df %>%
        filter(hugo_symbol %in% !!gene_set) %>%
        pull(gene_effect_scaled)
    y <- df %>%
        filter(!hugo_symbol %in% !!gene_set) %>%
        pull(gene_effect_scaled)

    if (length(gene_set) < 15 | length(gene_set) > 500) { return(NULL) }
    if (length(x) < 10 | length(y) < 20) { return(NULL) }

    wcrs_results <- tidy(wilcox.test(x, y, paried = FALSE))
    return(wcrs_results)
}

#### ---- Run the test on each cell line ---- ####

cellline_data <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    mutate(gene_effect_scaled = scale(gene_effect)[, 1]) %>%
    group_by(dep_map_id) %>%
    nest() %>%
    ungroup()

cell_line_info <- model_data %>%
    select(dep_map_id, cancer, allele) %>%
    unique()

cache("cellline_enrichment_results",
      depends = c("gene_sets_df", "cellline_data"),
{
    cellline_enrichment_results <- purrr::map(
        1:nrow(gene_sets_df),
        function(i) {
            gene_set_name <- unlist(gene_sets_df$gene_set[i])
            genes <- unlist(gene_sets_df$genes[i])

            results <- cellline_data %>%
                mutate(
                    gene_set = !!gene_set_name,
                    enrichment_results = purrr::map(data,
                                                    ranksum_enrichment,
                                                    gene_set = genes)
                ) %>%
                select(-data) %>%
                unnest(enrichment_results)
            return(results)
        }
    ) %>%
        bind_rows() %>%
        janitor::clean_names() %>%
        group_by(dep_map_id) %>%
        mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
        ungroup() %>%
        left_join(cell_line_info, by = "dep_map_id")
    return(cellline_enrichment_results)
})


#### ---- Run the test on cell lines of an allele ---- ####

allele_data <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol, allele) %>%
    summarise(avg_gene_effect = mean(gene_effect, rm.na = TRUE)) %>%
    group_by(cancer, hugo_symbol) %>%
    mutate(gene_effect_scaled = scale(avg_gene_effect)[, 1]) %>%
    group_by(cancer, allele) %>%
    nest() %>%
    ungroup()

cache("allele_enrichment_results",
      depends = c("gene_sets_df", "allele_data"),
{
    allele_enrichment_results <- purrr::map(
        1:nrow(gene_sets_df),
        function(i) {
            gene_set_name <- unlist(gene_sets_df$gene_set[i])
            genes <- unlist(gene_sets_df$genes[i])

            results <- allele_data %>%
                mutate(
                    gene_set = !!gene_set_name,
                    enrichment_results = purrr::map(data,
                                                    ranksum_enrichment,
                                                    gene_set = genes)
                ) %>%
                select(-data) %>%
                unnest(enrichment_results)
            return(results)
        }
    ) %>%
        bind_rows() %>%
        janitor::clean_names() %>%
        group_by(cancer, gene_set) %>%
        mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
        ungroup()
})


#### ---- Plotting allele results ---- ####

p <- allele_enrichment_results %>%
    filter(adjusted_p_value < 0.1) %>%
    mutate(norm_statistic = scales::rescale(statistic, to = c(-1, 1))) %>%
    ggplot(aes(x = allele, y = gene_set)) +
    facet_wrap(~ cancer, nrow = 1, scale = "free") +
    geom_point(
        aes(size = -log10(adjusted_p_value),
            color = norm_statistic
    )) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw(base_family = "arial", base_size = 7) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1.0),
        axis.title.x = element_blank()
    ) +
    labs(
        title = "Gene set enrichment in DepMap",
        y = "gene set",
        color = "test stat. (norm.)",
        size = "-log10( FDR-adj. p-value )"
    )
save_path <- plot_path("10_50_single-sample-rank-sum",
                       "enriched_gene_sets.svg")
ggsave_wrapper(p, save_path, "wide")


gene_set_waterfall_plot <- function(cancer, allele, gene_set,
                                    p_value, statistic, adjusted_p_value,
                                    ...) {
    genes <- gene_sets_df %>%
        filter(gene_set == !!gene_set) %>%
        pull(genes) %>%
        unlist()

    p_data <- allele_data %>%
        filter(cancer == !!cancer) %>%
        unnest(data) %>%
        filter(hugo_symbol %in% !!genes) %>%
        mutate(
            allele_gene = paste0(allele, "_", hugo_symbol),
            allele_gene = fct_reorder(allele_gene, gene_effect_scaled),
            y = gene_effect_scaled,
            ymin = purrr::map_dbl(y, ~ min(.x, 0)),
            ymax = purrr::map_dbl(y, ~ max(.x, 0)),
            label = ifelse(allele == !!allele, hugo_symbol, NA),
            label_hjust = ifelse(y > 0, 1.0, 0.0)
        )

    p_waterfall <- p_data %>%
        ggplot(aes(x = allele_gene)) +
        geom_linerange(
            aes(ymax = ymax, ymin = ymin),
            group = "allele_gene",
            size = 0.5,
            color = "grey50"
        ) +
        geom_point(
            aes(y = y, color = allele),
            size = 0.8
        ) +
        geom_text(
            aes(label = label, hjust = label_hjust),
            y = 0,
            angle = 60,
            family = "Arial",
            size = 1.5
        ) +
        scale_color_manual(values = short_allele_pal) +
        theme_bw(
            base_family = "Arial",
            base_size = 7
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = c(0.8, 0.2)
        ) +
        labs(
            title = glue("{cancer} KRAS {allele} enriched in {gene_set}"),
            y = "average scaled depletion effect"
        )

    p_tile <- p_data %>%
        ggplot(aes(x = allele_gene)) +
        geom_tile(aes(fill = allele, y = "KRAS"), color = NA) +
        scale_fill_manual(values = short_allele_pal) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_bw(
            base_family = "Arial",
            base_size = 7
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.position = "none"
        )

    p_arranged <- (p_waterfall / p_tile) +
        plot_layout(heights = c(40, 1))

    save_path <- plot_path(
        "10_50_single-sample-rank-sum",
        glue("waterfall_{cancer}_{allele}_{gene_set}.svg")
    )

    save_size <- case_when(
        length(gene_set) < 20 ~ c(5, 4),
        length(gene_set) < 50 ~ c(8, 4),
        TRUE ~ c(12, 4),
    )
    ggsave_wrapper(p_arranged, save_path,
                   width = save_size[[1]],
                   height = save_size[[2]])
}

allele_enrichment_results %>%
    filter(adjusted_p_value < 0.1) %>%
    pwalk(gene_set_waterfall_plot)
