# Some of the genes have disagreeing comutation interactions with multiple
# KRAS alleles - i.e. a single gene will have an increased comutation
# interaction with one allele, but a reduced comutation interaction with
# another allele. These will be referenced as "conflicts" below.
# Here I made bar plots of the comutation frequencies for each gene that has
# this behavior.


GRAPHS_DIR <- "20_41_disagreeing-interactions_mutfreq-barplot"
reset_graph_directory(GRAPHS_DIR)

# A gene in conflict is one with an interaction with at least two alleles
# and has multiple types of interaction.
conflicts_df <- genetic_interaction_df %>%
    group_by(hugo_symbol, cancer) %>%
    filter(n_distinct(allele) > 1 & n_distinct(genetic_interaction) == 2) %>%
    ungroup()



comut_summary_df <- fisher_comut_df %>%
    mutate(allele = str_remove_all(kras_allele, "KRAS_")) %>%
    select(cancer, allele, hugo_symbol, odds_ratio, n00, n10, n01, n11)


make_mutfreq_barplots <- function(cancer, data, ...) {
    p <- data %>%
        mutate(
            log_or = log10(odds_ratio),
            allele = factor_alleles(allele),
            genetic_interaction = switch_comut_terms(genetic_interaction)
        ) %>%
        ggplot(aes(x = allele, y = log_or)) +
        facet_wrap(~ hugo_symbol, scales = "free_x") +
        geom_col(aes(fill = genetic_interaction)) +
        geom_hline(yintercept = 0, size = 0.5, color = "black") +
        scale_fill_manual(values = comut_updown_pal) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            axis.title.x = element_blank(),
            strip.background = element_blank(),
            legend.position = "bottom"
        )  +
        labs(
            y = "log(OR)",
            fill = "comutation"
        )
    return(p)
}


save_mutfreq_barplots <- function(cancer, data, bar_plot, ...) {
    ggsave_wrapper(
        bar_plot,
        plot_path(GRAPHS_DIR, glue("mutfreq-barplot_{cancer}.svg")),
        "medium"
    )
}



conflicts_df %>%
    select(cancer, allele, hugo_symbol, genetic_interaction, p_val) %>%
    left_join(comut_summary_df, by = c("cancer", "allele", "hugo_symbol")) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(bar_plot = purrr::map2(cancer, data, make_mutfreq_barplots)) %>%
    pmap(save_mutfreq_barplots)
