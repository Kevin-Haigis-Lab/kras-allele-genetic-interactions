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
        legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6)
    ) +
    labs(tag = "a")


#### ---- B. Distribution of KRAS alleles ---- ####
# The distribution of KRAS alleles across cancers and codon.
# original script: "src/##_##_ORIGINAL-SCRIPT.R"

panel_B <- plot_spacer()


#### ---- C. Percent of samples with KRAS mutation ---- ####
# Percent of samples with a KRAS mutation per cancer.
# original script: "src/##_##_ORIGINAL-SCRIPT.R"

panel_C <- plot_spacer()


#### ---- D. Mutational signature demonstration ---- ####
# Demonstration of mutational signatures decomposing a mutational spectrum.
# original script: "src/##_##_ORIGINAL-SCRIPT.R"

panel_D <- plot_spacer()


#### ---- E. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/##_##_ORIGINAL-SCRIPT.R"

panel_E <- plot_spacer()


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    row_1 <- panel_A  | panel_B | panel_C
        plot_layout(widths = c(1, 3, 1))

    full_figure <- row_1 / panel_D / panel_E

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
