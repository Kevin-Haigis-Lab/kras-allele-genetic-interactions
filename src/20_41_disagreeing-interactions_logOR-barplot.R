# Some of the genes have disagreeing comutation interactions with multiple
# KRAS alleles - i.e. a single gene will have an increased comutation
# interaction with one allele, but a reduced comutation interaction with
# another allele. These will be referenced as "conflicts" below.
# Here I made bar plots of the comutation frequencies for each gene that has
# this behavior.


GRAPHS_DIR <- "20_41_disagreeing-interactions_logOR-barplot"
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


make_logOR_barplots <- function(cancer, data, ...) {
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


save_logOR_barplots <- function(cancer, data, bar_plot, ...) {
    ggsave_wrapper(
        bar_plot,
        plot_path(GRAPHS_DIR, glue("logOR-barplot_{cancer}.svg")),
        "medium"
    )
}



conflicts_df %>%
    select(cancer, allele, hugo_symbol, genetic_interaction, p_val) %>%
    left_join(comut_summary_df, by = c("cancer", "allele", "hugo_symbol")) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(bar_plot = purrr::map2(cancer, data, make_logOR_barplots)) %>%
    pmap(save_logOR_barplots)



#### ---- Saving specific interactoins for figures ---- ####

# A S4 class to hold information for saving specific interactions for figures.
setClass("SpecificInteraction",
    slots = c(
        cancer = "character",
        genes =  "character",
        fig_num = "numeric",
        supp = "logical"
    ),
    prototype = list(
        cancer = NA_character_,
        genes = NA_character_,
        fig_num = NA_real_,
        supp = FALSE
    )
)


# A list of interactions to save.
specific_interactions <- c(
    new("SpecificInteraction",
        cancer = "PAAD",
        genes = c("TP53", "RNF43", "MAP2K4", "RBM10"),
        fig_num = 12,
        supp = TRUE
    )
)


# Save the ggplot object for a figure.
save_logOR_ggproto <- function(bar_plot, si_obj, ...) {
    saveRDS(
        bar_plot,
        get_fig_proto_path(
            glue("log-odds-ratio_barplot_{si_obj@cancer}"),
            figure_num = si_obj@fig_num,
            supp = si_obj@supp
        )
    )

}


for (si in specific_interactions) {
    conflicts_df %>%
        select(cancer, allele, hugo_symbol, genetic_interaction, p_val) %>%
        filter(cancer == si@cancer & hugo_symbol %in% si@genes) %>%
        left_join(
            comut_summary_df,
            by = c("cancer", "allele", "hugo_symbol")
        ) %>%
        mutate(hugo_symbol = factor(hugo_symbol, levels = si@genes)) %>%
        group_by(cancer) %>%
        nest() %>%
        mutate(bar_plot = purrr::map2(cancer, data, make_logOR_barplots)) %>%
        pmap(save_logOR_ggproto, si_obj = si)
}
