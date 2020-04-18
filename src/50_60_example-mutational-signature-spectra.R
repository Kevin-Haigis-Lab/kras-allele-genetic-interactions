# An example of the mutational spectra and decomposition into signatures.

GRAPHS_DIR <- "50_60_example-mutational-signature-spectra"
reset_graph_directory(GRAPHS_DIR)


coad_genome_muts <- trinucleotide_mutations_df %>%
    filter(!is_hypermutant & cancer == "COAD" & target == "genome")

style_mutational_spectrum_plot <- function(p) {
    p +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = mutsig_context_group_pal) +
        # scale_fill_brewer(type = "seq", palette = "PuBuGn") +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
            axis.title.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.position = "none",
            panel.grid.major.x = element_blank()
        )
}


extract_context_group <- function(x) {
    str_extract(x, "(?<=\\[).*(?=\\])")
}


#### ---- Mutational spectrum of COAD WHS samples ---- ####

coad_mut_spectrum <- coad_genome_muts %>%
    count(context, tricontext) %>%
    mutate(context_grp = extract_context_group(tricontext),
           mutation = str_remove(tricontext, "\\[[CT]\\>"),
           mutation = str_remove(mutation, "\\]"),
           n_norm = n / sum(n)) %>%
    ggplot(aes(x = context, y = n_norm)) +
    facet_grid(. ~ context_grp, scales = "free_x") +
    geom_col(aes(fill = context_grp)) +
    labs(y = "distribution",
         title = "COAD mutational spectrum")
coad_mut_spectrum <- style_mutational_spectrum_plot(coad_mut_spectrum)
ggsave_wrapper(
    coad_mut_spectrum,
    plot_path(GRAPHS_DIR, "coad_mut_spectrum.svg"),
    "wide"
)
saveFigRds(coad_mut_spectrum, "coad_mut_spectrum")



#### ---- Example signature spectra ---- ####

example_signatures <- c("1", "4", "5", "8", "9", "18")

example_mutsig_spectra <- mutational_signature_spectra %>%
    filter(signature %in% example_signatures) %>%
    mutate(signature = paste("sig.", signature),
           signature = fct_inorder(signature),
           context = str_remove(tricontext, "\\>[:alpha:]\\]"),
           context = str_remove(context, "\\["),
           context_grp = extract_context_group(tricontext)) %>%
    ggplot(aes(x = context, y = composition)) +
    facet_grid(signature ~ context_grp, scales = "free") +
    geom_col(aes(fill = context_grp)) +
    labs(y = "composition",
         title = "Example signature mutational spectra")
example_mutsig_spectra <- style_mutational_spectrum_plot(example_mutsig_spectra)
ggsave_wrapper(
    example_mutsig_spectra,
    plot_path(GRAPHS_DIR, "example_mutsig_spectra.svg"),
    "wide"
)
saveFigRds(example_mutsig_spectra, "example_mutsig_spectra")


# Also save each signature as a separate object.
for (sig in example_signatures) {
    eg_spectrum <- mutational_signature_spectra %>%
        filter(signature %in% !!sig) %>%
        mutate(signature = paste("sig.", signature),
               signature = fct_inorder(signature),
               context = str_remove(tricontext, "\\>[:alpha:]\\]"),
               context = str_remove(context, "\\["),
               context_grp = extract_context_group(tricontext)) %>%
        ggplot(aes(x = context, y = composition)) +
        facet_grid(. ~ context_grp, scales = "free") +
        geom_col(aes(fill = context_grp)) +
        labs(y = "composition",
             title = "Example signature mutational spectra")
    saveFigRds(style_mutational_spectrum_plot(eg_spectrum),
               as.character(glue("sig{sig}_example_mutsig_spectra")))
}
