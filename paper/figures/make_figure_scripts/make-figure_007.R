# Figure 007. Distribution of mutational signatures across alleles

FIGNUM <- 7

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig7 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm"))
    )
}


#### ---- A. Distirubiton of mutational signatures in each sample ---- ####
# The distribution of mutational signature levels in each sample (barplot).
# original script: "src/50_20_mutsignatures-distributions.R"

panel_A <- read_fig_proto("signature-level-per-sample") +
    scale_fill_manual(
        values = mutsig_descrpt_pal,
        guide = guide_legend(
            nrow = 1,
            label.position = "top",
            label.hjust = 0.5,
            label.vjust = -11
        )
    ) +
    labs(
        y = "signature level",
        fill = "signature",
        tag = "a"
    ) +
    theme_fig7() +
    theme(
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.25),
        axis.text.x = element_blank(),
        legend.key.size = unit(4, "mm")
    )


#### ---- B. Distirubiton of mutational signature levels ---- ####
# The distribution of mutational signature levels in each sample (boxplot).
# original script: "src/50_20_mutsignatures-distributions.R"

panel_B <- read_fig_proto("signature-level-boxplots_with0") +
    scale_y_continuous(expand = expansion(add = c(0.02, 0.02))) +
    labs(
        tag = "b",
        y = "signature level"
    ) +
    theme_fig7() +
    theme(
        legend.position = "none"
    )


#### ---- C. Distirubiton of mutational signatures by allele ---- ####
# The distribution of mutational signature levels in samples with each allele.
# original script: "src/50_20_mutsignatures-distributions.R"

panel_C <- read_fig_proto("mutational-signatures-distribution-by-allele") +
    labs(tag = "c") +
    theme_fig7() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank()
    )


#### ---- D. Levels of clock vs. non-clock mutational signatures ---- ####
# The levels of clock and non-clock mutational signatures per cancer.
# original script: "src/50_20_mutsignatures-distributions.R"

panel_D <- read_fig_proto("clock-signatures_violin-box") +
    labs(tag = "d") +
    theme_fig7() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank()
    )


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (panel_A + panel_B + plot_layout(widths = c(1))) /
        (panel_C + panel_D + plot_layout(widths = c(7, 3))) /
        guide_area() +
        plot_layout(heights = c(15, 10, 1), guides = "collect")

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
