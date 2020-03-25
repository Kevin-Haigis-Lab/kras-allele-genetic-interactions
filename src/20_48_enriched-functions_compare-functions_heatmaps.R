# Make heatmaps of genes in multiple gene sets.

GRAPHS_DIR <- "20_48_enriched-functions_compare-functions_heatmaps"
reset_graph_directory(GRAPHS_DIR)


#### ---- Input ---- ####
# Add to the data frame below to make new plots.
# Each group `grp` will be plotted, separately.
geneset_groups <- tibble::tribble(
    ~grp, ~cancer, ~datasource, ~term,
    1, "PAAD", "Transcription_Factor_PPIs", "SMAD1",
    1, "PAAD", "Transcription_Factor_PPIs", "SMAD2",
    1, "PAAD", "Transcription_Factor_PPIs", "SMAD3",
)


#### ---- Subroutines ---- ####

# Extract data from the Enrichr results by passing a data frame with the
# columns to use for a right join with Enrichr.
get_enrichr_results <- function(df, p_value_max = 0.05, min_overlap = 3) {
    # Check that the necessary column names are present.
    df %>% verify(has_all_names("cancer", "term", "datasource"))

    # Extract the relevant results from `enrichr_tib`.
    enrichr_tib %>%
        select(-gene_list) %>%
        unnest(enrichr_res) %>%
        right_join(df, by = c("cancer", "datasource", "term")) %>%
        mutate(n_overlap = get_enrichr_overlap_int(overlap)) %>%
        filter(
            adjusted_p_value < !!p_value_max & n_overlap >= !!min_overlap
        ) %>%
        mutate(overlap_genes = str_split(genes, ";")) %>%
        select(-genes) %>%
        mutate(
            mod_term = standardize_enricher_terms(term),
            mod_term = map2_chr(term, datasource, mod_term_for_datasource)
        )
}

# A data frame of summary comutation data.
comut_summary_df <- fisher_comut_df %>%
    mutate(allele = str_remove_all(kras_allele, "KRAS_")) %>%
    select(cancer, allele, hugo_symbol, odds_ratio, n00, n10, n01, n11)

# Get comutation data. Leave any parameter as `NULL` to not filter.
get_comutation_results <- function(cancers = NULL,
                                   alleles = NULL,
                                   genes = NULL) {
    comut_summary_df %>%
        filter(
            (cancer %in% !!cancers | is.null(!!cancers)) &
            (allele %in% !!alleles | is.null(!!alleles)) &
            (hugo_symbol %in% !!genes | is.null(!!genes)))
}

# Make the heatmap for whether the gene is in a geneset or not.
make_binary_geneset_heatmap <- function(df) {
    p <- df %>%
        mutate(a = TRUE) %>%
        tidyr::complete(nesting(datasource, term), hugo_symbol,
                        fill = list(a = FALSE)) %>%
        mutate(term = factor(term, levels = sort(unique(term)))) %>%
        ggplot(aes(x = hugo_symbol, y = term)) +
        geom_tile(aes(fill = a), color = "grey25") +
        scale_fill_manual(values = c("white", "black")) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none",
            axis.ticks = element_blank(),
            panel.grid = element_blank()
        )
    return(p)
}

# Sort genes by their comutation frequency and allele. ("memosort")
sort_genes_by_comutfreq <- function(df) {
    allele_order <- levels(df$allele)
    allele_order <- allele_order[allele_order %in% unique(df$allele)]

    gene_order <- df %>%
        mutate(
            allele_score = map_dbl(allele, ~ which(allele_order == .x)),
            score = (allele_score**3) * comut_freq
        ) %>%
        group_by(hugo_symbol) %>%
        summarise(score = sum(score)) %>%
        ungroup() %>%
        arrange(-score) %>%
        pull(hugo_symbol)

    df %>%
        mutate(hugo_symbol = factor(hugo_symbol, levels = gene_order))
}

# Make the heatmap of log(OR) or stacked bar-plot showing comutation rates.
make_comutation_heatmap <- function(df,
                                    which_metric = c("OR", "comut_freq"),
                                    max_val = 1) {
    which_metric <- which_metric[[1]]
    df <- df %>%
        select(-datasource, -term) %>%
        unique()
    if (which_metric == "OR") {
        p <- df %>%
            mutate(
                log_or = log10(odds_ratio),
                log_or = scales::squish(log_or, range = c(-max_val, max_val))
            ) %>%
            ggplot(aes(x = hugo_symbol, y = allele)) +
            geom_tile(aes(fill = log_or), color = "grey50") +
            scale_fill_gradient2(
                high = comut_updown_pal["increased"],
                low = comut_updown_pal["reduced"],
                mid = "white"
            ) +
            scale_y_discrete(expand = c(0, 0))
    } else if (which_metric == "comut_freq") {
        mod_df <- df %>%
            mutate(comut_freq = n11 / (n11 + n01)) %>%
            group_by(hugo_symbol) %>%
            mutate(comut_freq = comut_freq / sum(comut_freq)) %>%
            ungroup() %>%
            sort_genes_by_comutfreq()
        p <- mod_df %>%
            ggplot(aes(x = hugo_symbol, y = comut_freq)) +
            geom_col(aes(fill = allele), position = "fill") +
            scale_fill_manual(values = short_allele_pal) +
            scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.0, 0.25))
    } else {
        stop(as.character(glue("Do not recognize the metric '{which_metric}'")))
    }

    p <- p +
        scale_x_discrete(expand = c(0, 0)) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.key.size = unit(2, "mm")
        ) +
        labs(
            fill = "log10(OR)",
            y = "distribution of\ncomutation events"
        )
    return(p)
}

# Make the heatmap / stacked bar-plot for a set of terms and datasources.
make_comparison_heatmap <- function(grp, data, ...) {
    enrichr_data_df <- get_enrichr_results(data)
    geneset_genes_df <- enrichr_data_df %>%
        select(datasource, term, overlap_genes) %>%
        unnest(overlap_genes) %>%
        unique() %>%
        dplyr::rename(hugo_symbol = overlap_genes)

    all_genes <- unique(unlist(enrichr_data_df$overlap_genes))
    comut_df <- get_comutation_results(
        cancers = unique(data$cancer),
        alleles = unique(enrichr_data_df$allele),
        genes = all_genes
    )

    all_data <- left_join(comut_df, geneset_genes_df, by = "hugo_symbol") %>%
        mutate(allele = factor_alleles(allele))

    top_hm <- make_binary_geneset_heatmap(all_data)
    bottom_hm <- make_comutation_heatmap(all_data, which_metric = "comut_freq")

    hm <- (top_hm / bottom_hm) + plot_layout(heights = c(1, 3))
    return(hm)
}

# Save info for a figure.
save_fig_info <- tibble::tribble(
    ~grp, ~info,
       1,  list(fig_num = 12, supp = TRUE),
)

# Save the comparison heatmaps.
# Add the information to `save_fig_info` to save for a figure.
save_comparison_heatmaps <- function(grp, data, hm, ...) {
    cancer <- unique(data$cancer)
    save_name <- as.character(glue("comparison-heatmap_{cancer}-{grp}.svg"))
    ggsave_wrapper(
        hm,
        plot_path(GRAPHS_DIR, save_name),
        width = 4, height = 2
    )

    # Save to figure folder.
    if (grp %in% save_fig_info$grp) {
        info <- save_fig_info %>%
            filter(grp == !!grp) %>%
            pull(info) %>%
            unlist(recursive = FALSE)
        saveFigRds(hm, save_name)
    }
}

# Make the comparison heatmaps / stacked bar-plots for the selected groups
# of gene sets.
geneset_groups %>%
    group_by(grp) %>%
    nest() %>%
    mutate(hm = purrr::map2(grp, data, make_comparison_heatmap)) %>%
    pwalk(save_comparison_heatmaps)
