# What is the probability that each mutational signature caused the allele
# of each sample.

library(patchwork)

GRAPHS_DIR <- "50_30_mutsignatures_prob-causing-allele"
reset_graph_directory(GRAPHS_DIR)


# A data frame of the number of samples with each allele.
alleles_frequency_per_cancer_df <- mutational_signatures_df %>%
    mutate(ras_allele = str_remove_all(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, cancer, ras_allele) %>%
    unique() %>%
    count(cancer, ras_allele)

# Which alleles to plot for a cancer based off of the number of samples
# available with the allele.
alleles_to_plot <- function(cancer, min_num = 15) {
    alleles_frequency_per_cancer_df %>%
        filter(cancer == !!cancer & n >= min_num) %>%
        pull(ras_allele) %>%
        unique()
}


# Just the columns with mutational signature information.
mutsig_noartifact_df_select <- mutsig_noartifact_df %>%
    select(tumor_sample_barcode, signature, contribution, description)


# A compilation of the mutational signature information for the KRAS mutations.
#   contribution: the level of the signature in the sample
#   composition: the level of the context in the signature
#   casuation_prob: the likelihood of the mutation coming from the
#                   signature's mutagen
kras_allele_mutsig_df <- trinucleotide_mutations_df %>%
    filter(hugo_symbol == "KRAS" &
           kras_allele != "WT" &
           cancer != "SKCM" &
           !is_hypermutant) %>%
    left_join(mutational_signature_spectra, by = "tricontext") %>%
    left_join(mutsig_noartifact_df_select,
              by = c("tumor_sample_barcode", "signature")) %>%
    filter(!is.na(description)) %>%
    mutate(causation_prob = composition * contribution)


# Creates a column bar plot with position "fill" of the probability that
# the KRAS mutation was from a specific mutational signature.
plot_probability_of_causation <- function(cancer, data,
                                          min_allele_num = 15, ...) {
    # Alleles to include in the plots.
    alleles <- alleles_to_plot(cancer, min_allele_num)

    # Prepare the data for plotting.
    plot_data <- data %>%
        filter(kras_allele %in% !!alleles) %>%
        group_by(kras_allele, description) %>%
        summarise(
            prob_sum = sum(causation_prob),
            prob_avg = mean(causation_prob)
        ) %>%
        group_by(kras_allele) %>%
        mutate(prob_sum = prob_sum / sum(prob_sum)) %>%
        ungroup()

    # Plot.
    p <- plot_data %>%
        mutate(
            kras_allele = factor(kras_allele, names(short_allele_pal)),
            description = factor(description, levels = names(mutsig_descrpt_pal))
        ) %>%
        ggplot(aes(x = kras_allele, y = prob_sum)) +
        geom_col(aes(fill = description), position = "fill") +
        scale_fill_manual(
            values = mutsig_descrpt_pal,
            guide = guide_legend(
                nrow = 1,
                title.position = "left",
                label.position = "top",
                label.hjust = 0.5,
                label.vjust = -5
            )
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw(
            base_size = 9,
            base_family = "Arial"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "bottom",
            legend.key.size = unit(4, "mm"),
            legend.spacing.y = unit(0, "mm")
        ) +
        labs(
            title = cancer,
            y = "prob. of causing allele"
        )
    return(p)
}

barplots <- kras_allele_mutsig_df %>%
    group_by(cancer) %>%
    nest() %>%
    arrange(cancer) %>%
    pmap(plot_probability_of_causation)
combined_plots <-  wrap_plots(barplots) / guide_area() +
    plot_layout(guides = 'collect', heights = c(14, 1))
ggsave_wrapper(
    combined_plots,
    plot_path(GRAPHS_DIR, "probability-mutsig-caused-allele.svg"),
    'wide'
)
