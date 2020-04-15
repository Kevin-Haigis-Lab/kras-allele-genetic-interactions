# Figure 028. #> Mutational signature main figure.
# (Based off of sketch from PJP.)

FIGNUM <- 28

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig28 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. KRAS hot-spot mutation frequency barplot ---- ####
# The distribution of mutations at the 4 hotspots.
# original script: "src/90_05_kras-allele-distribution.R"

panel_A <- read_fig_proto("lollipop-kras_hotspot-only") +
    theme_fig28() +
    theme(
        panel.grid.major.x = element_blank(),
        axis.title.y = element_markdown(),
        legend.position = c(0.80, 0.82),
        legend.background = element_blank(),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6),
        plot.margin = margin(0, 2, 0, 0, "mm")
    ) +
    labs(tag = "a",
         y = "number of samples (*log*<sub>10</sub>)")


#### ---- B. Distribution of KRAS alleles ---- ####
# The distribution of KRAS alleles across cancers and codon.
# original script: "src/90_05_kras-allele-distribution.R"

panel_B <- read_fig_proto("allele_dist_dotplot") +
    theme_fig28() +
    theme(
        axis.ticks = element_blank(),
        axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold")
    ) +
    labs(tag = "b")


#### ---- C. Percent of samples with KRAS mutation ---- ####
# Percent of samples with a KRAS mutation per cancer.
# original script: "src/90_05_kras-allele-distribution.R"

panel_C <- read_fig_proto("cancer_freq_kras_mut_column") +
    scale_x_continuous(
        expand = expansion(add = c(0, 0.02)),
        breaks = c(0.2, 0.4, 0.6, 0.8),
        labels = function(x) { paste0(round(x * 100), "%") }
    ) +
    theme_fig28() +
    theme(
        axis.title.x = element_markdown(),
        axis.text.x = element_text(angle = -90),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank()
    ) +
    labs(tag = "c")


#### ---- D. Mutational signature demonstration ---- ####
# Demonstration of mutational signatures decomposing a mutational spectrum.
# original script: "src/50_60_example-mutational-signature-spectra.R"

mutational_spectrum_theme <- function(...) {
    theme_fig28(...) %+replace%
    theme(
        plot.title = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(size = 6, angle = -90,
                                   hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        legend.position = "none",
        panel.grid.major.x = element_blank()
    )
}

panel_D_1 <- read_fig_proto("coad_mut_spectrum") +
    mutational_spectrum_theme() +
    labs(title = "COAD mutational spectrum")

panel_D_2 <- read_fig_proto("example_mutsig_spectra") +
    mutational_spectrum_theme() +
    labs(title = NULL)

panel_D <- (
        (plot_spacer() / panel_D_1 / plot_spacer()) +
        plot_layout(heights = c(1, 2, 1))
    ) |
    panel_D_2


#### ---- E. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/##_##_ORIGINAL-SCRIPT.R"

panel_E <- plot_spacer()


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    row_1 <- (panel_A  | panel_B | panel_C) +
        plot_layout(widths = c(2, 5, 1))

    full_figure <- (row_1 / panel_D / panel_E) +
        plot_layout(heights = c(2, 3, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
