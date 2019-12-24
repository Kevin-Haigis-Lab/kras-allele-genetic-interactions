# What is the probability that each mutational signature caused the allele
# of each sample.

GRAPHS_DIR <- "50_30_mutsignatures_prob-causing-allele"
reset_graph_directory(GRAPHS_DIR)


# A data frame of the number of samples with each allele.
alleles_frequency_per_cancer_df <- mutational_signatures_df %>%
    filter(!is_hypermutant) %>%
    mutate(kras_allele = str_remove_all(ras_allele, "KRAS_")) %>%
    select(tumor_sample_barcode, cancer, kras_allele) %>%
    unique() %>%
    count(cancer, kras_allele)

# Which alleles to plot for a cancer based off of the number of samples
# available with the allele.
alleles_to_plot <- function(cancer, min_num = 15) {
    alleles_frequency_per_cancer_df %>%
        filter(cancer == !!cancer & n >= min_num) %>%
        pull(kras_allele) %>%
        unique()
}


# Just the columns with mutational signature information.
mutsig_noartifact_df_select <- mutsig_noartifact_df %>%
    filter(!is_hypermutant) %>%
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
            kras_allele = factor_alleles(kras_allele),
            description = factor(description, levels = names(mutsig_descrpt_pal))
        ) %>%
        ggplot(aes(x = kras_allele, y = prob_sum)) +
        geom_col(aes(fill = description), position = "fill") +
        scale_fill_manual(
            values = mutsig_descrpt_pal,
            guide = guide_legend(
                nrow = 1,
                title.position = "left",
                title.vjust = 0.2,
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
            y = "probability",
            fill = "signature"
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
saveRDS(
    combined_plots,
    get_fig_proto_path("probability-mutsig-caused-allele", 1)
)

#### ---- Distribution of select signatures ---- ####


plot_signature_probability <- function(cancer, signature, min_allele_num = 15) {
    # Alleles to include in the plots.
    alleles <- alleles_to_plot(cancer, min_allele_num)

    # Prepare the data for plotting.
    plot_data <- kras_allele_mutsig_df %>%
        filter(
            !is_hypermutant &
            cancer == !!cancer &
            signature %in% !!signature &
            kras_allele %in% !!alleles
        )

    p <- plot_data %>%
        mutate(
            kras_allele = factor_alleles(kras_allele, reverse = TRUE)
        ) %>%
        ggplot(aes(x = kras_allele, y = causation_prob)) +
        ggbeeswarm::geom_quasirandom(
            aes(color = kras_allele),
            alpha = 0.7,
            size = 0.15,
            varwidth = FALSE,
            method = "quasirandom"
        ) +
        scale_color_manual(
            values = short_allele_pal,
            guide = FALSE
        ) +
        scale_fill_manual(
            values = short_allele_pal,
            guide = FALSE
        ) +
        scale_y_continuous(
            expand = expand_scale(mult = c(0.01, 0.05))
        ) +
        theme_bw(
            base_size = 8,
            base_family = "Arial"
        ) +
        theme(
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = glue("Sig. {signature} in {cancer}"),
            y = "signature level"
        )

    return(p)
}

select_sig_plots <- tibble::tribble(
    ~cancer, ~signature,
     "COAD",       "18",
     "LUAD",        "4",
       "MM",        "9",
     "PAAD",        "8"
) %>%
    pmap(plot_signature_probability)
select_sig_combined_plots <-  wrap_plots(select_sig_plots)
ggsave_wrapper(
    select_sig_combined_plots,
    plot_path(GRAPHS_DIR, "contribution-of-select-signatures.svg"),
    width = 5, height = 4
)

saveRDS(
    select_sig_combined_plots,
    get_fig_proto_path("contribution-of-select-signatures", 1)
)
