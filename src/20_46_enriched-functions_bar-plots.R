# Another layer of abstraction on top of the functional enrichment results.
# Bar-plots showing the rate of comutation of all the genes in the function
# in one allele compared to the rest samples.


GRAPHS_DIR <- "20_46_enriched-functions_bar-plots"
reset_graph_directory(GRAPHS_DIR)

# Get all of the genes for a data source and term used by Enrichr.
# Make sure to use the original data source name and term.
genes_for <- function(datasource, term) {
    genes <- enrichr_genesets %>%
        filter(datasource == !!datasource & term == !!term) %>%
        u_pull(genes)
    return(unlist(genes))
}


# A data frame with the number of samples per allele per cancer.
num_samples_per_allele <- cancer_coding_muts_df %>%
    filter(!is_hypermutant & cancer != "SKCM") %>%
    group_by(cancer, ras_allele) %>%
    summarise(num_cancer_samples = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    mutate(allele = str_remove_all(ras_allele, "KRAS_"))

# Get the number of samples of a given allele and not in `alleles_to_remove`.
# If `allele == "Other"`: then only `alleles_to_remove` are removed.
get_number_of_samples <- function(cancer, allele, not_allele = FALSE) {
    muts_tib <- num_samples_per_allele %>%
        filter(cancer %in% !!cancer)

    if (not_allele) {
        muts_tib %<>% filter(allele != !!allele)
    } else  {
        muts_tib %<>% filter(allele == !!allele)
    }

    sum(unlist(muts_tib$num_cancer_samples))
}
get_number_of_samples <- memoise::memoise(get_number_of_samples)


# Calculate the mutational frequency of `genes` in the samples with an `allele`
#   in `cancer` but not with the allele in `alleles_to_remove`.
# If `allele == "Other"`: then only `alleles_to_remove` are removed.
calc_mut_freq <- function(cancer, allele, allele_grp, direction_genes, ...) {
    muts_tib <- cancer_coding_muts_df %>%
        filter(!is_hypermutant) %>%
        filter(
            cancer == !!cancer & hugo_symbol %in% !!unlist(direction_genes)
        ) %>%
        mutate(allele = str_remove_all(ras_allele, "KRAS_"))

    if (allele_grp == "Other") {
        muts_tib %<>% filter(allele != !!allele)
    } else {
        muts_tib %<>% filter(allele == !!allele)
    }

    n_mutants <- muts_tib %>%
        pull(tumor_sample_barcode) %>%
        unlist() %>%
        n_distinct()

    n_samples <- get_number_of_samples(cancer, allele, allele_grp == "Other")

    return(n_mutants / n_samples)
}



# Add columns and rows for allele group (the allele and "Other") and
#  the comutation direction ("increased" and "reduced").
add_new_group_columns_and_rows <- function(tib) {
    tib %>%
        mutate(allele_grp = purrr::map(allele, ~ c(.x, "Other"))) %>%
        unnest(allele_grp) %>%
        mutate(
            comut_direction = purrr::map(allele, ~ c("increased", "reduced"))
        ) %>%
        unnest(comut_direction)
}


# Filter the `overlap_genes` for whether the genes have increased
# or decreased rates of comutation.
genes_by_comut_direction <- function(cancer,
                                     allele,
                                     overlap_genes,
                                     comut_direction,
                                     ...) {
    # Get term depending on "increased" or "reduced" for comutation direction.
    comut_term <- ifelse(
        comut_direction == "increased", "comutation", "exclusivity"
    )

    # Get terms in `overlap_genes` that are in the correct comutation direction.
    genes <- genetic_interaction_df %>%
        filter(cancer == !!cancer & allele %in% !!allele) %>%
        filter(genetic_interaction == !!comut_term) %>%
        filter(hugo_symbol %in% !!overlap_genes) %>%
        u_pull(hugo_symbol)

    return(genes)
}


enrichr_comut_freq_tib <- enrichr_tib %>%
    select(-gene_list) %>%
    unnest(cols = enrichr_res) %>%
    mutate(n_overlap = get_enrichr_overlap_int(overlap)) %>%
    filter(adjusted_p_value < 0.05 & n_overlap >= 3) %>%
    mutate(overlap_genes = str_split(genes, ";")) %>%
    select(-genes) %>%
    mutate(term_mod = standardize_enricher_terms(term)) %>%
    select(cancer, allele, datasource, term, overlap_genes) %>%
    add_new_group_columns_and_rows()

enrichr_comut_freq_tib$direction_genes <- pmap(enrichr_comut_freq_tib,
                                               genes_by_comut_direction)
enrichr_comut_freq_tib$mut_freq <- pmap_dbl(enrichr_comut_freq_tib,
                                            calc_mut_freq)


# Order the terms be the difference in mutational frequency.
order_terms <- function(df) {
    term_order <- df %>%
        mutate(
            mut_freq = ifelse(allele_grp == "Other", -mut_freq, mut_freq)
        ) %>%
        group_by(term) %>%
        summarise(diff_total = sum(mut_freq)) %>%
        ungroup() %>%
        arrange(diff_total) %>%
        pull(term)
    df$term <- factor(df$term, levels = term_order)
    return(df)
}


# Bar plot of the enriched functions.
enriched_fxns_comutation_barplot <- function(data) {
    p <- data %>%
        mutate(
            mut_freq = ifelse(
                comut_direction == "reduced", -mut_freq, mut_freq
            ),
            term = str_wrap(term, 50),
            alpha_grp = ifelse(allele_grp == "Other", "other", "allele"),
            alpha_grp = factor(alpha_grp, levels = c("other", "allele"))
        ) %>%
        order_terms() %>%
        ggplot(aes(x = term, y = mut_freq)) +
        geom_col(
            aes(fill = comut_direction, alpha = alpha_grp),
            position = "stack"
        ) +
        geom_hline(
            yintercept = 0,
            size = 1,
            color = "grey25"
        ) +
        coord_flip() +
        scale_fill_manual(
            values = comut_updown_pal
        ) +
        scale_alpha_manual(
            values = c("other" = 0.4, "allele" = 0.95)
        ) +
        scale_y_continuous(
            label = abs
        ) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            axis.title.y = element_blank(),
            legend.position = "bottom",
            strip.background = element_blank()
        ) +
        labs(
            y = expression("reduced comutation" %<-% " " %->% "increased comutation"),
            alpha = "KRAS",
            fill = "comutation direction"
        )
    return(p)
}


# Save the plots to `GRAPHS_DIR`.
save_to_graphs_dir <- function(cancer,
                               allele = "AllAlleles",
                               datasource = "AllSources",
                               comut_plot,
                               data,
                               ...,
                               save_for_fig = NA,
                               supp = FALSE,
                               size = "large") {
    save_name <- as.character(
        glue("comut-barplot_{cancer}_{allele}_{datasource}.svg")
    )
    ggsave_wrapper(comut_plot, plot_path(GRAPHS_DIR, save_name), size)

    # Save the ggplot object for figure `save_for_fig`.
    if (!is.na(save_for_fig)) {
        saveFigRds(comut_plot, save_name)
    }
}


enrichr_comut_freq_tib %>%
    mutate(
        term = standardize_enricher_terms(term),
        term = map2_chr(term, datasource, mod_term_for_datasource)
    ) %>%
    group_by(cancer, allele, datasource) %>%
    nest() %>%
    mutate(
        comut_plot = purrr::map(data, enriched_fxns_comutation_barplot)
    ) %>%
    pwalk(save_to_graphs_dir)


# Read in a data frame (tibble) with a manually-selected subset of
# enriched functions.
source(file.path("src", "20_44_select-enriched-functions.R"))

additional_adjustments <- function(p) {
    p <- p +
        facet_grid(allele ~ .,
                   scales = "free_y",
                   space = "free_y")
    return(p)
}

enrichr_comut_freq_tib %>%
    mutate(term = standardize_enricher_terms(term)) %>%
    right_join(selected_enrichments,
               by = c("cancer", "datasource", "term")) %>%
    mutate(
        term = map2_chr(term, datasource, mod_term_for_datasource)
    ) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(
        comut_plot = purrr::map(data, enriched_fxns_comutation_barplot),
        comut_plot = purrr::map(comut_plot, additional_adjustments)
    ) %>%
    mutate(save_for_fig = ifelse(cancer == "LUAD", 3, NA)) %>%
    pwalk(save_to_graphs_dir, size = "medium")
