# Figure 008. Distribution of mutation frequency for genes that were
# removed due to low expression.

FIGNUM <- 8

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")
FIG_DIMENSIONS$height <- 100  # specific change for this figure


theme_fig8 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.position = "none",
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


#### ---- A. Distirubiton of mutational signatures in each sample ---- ####

# The distribution of frequency of mutations in genes removed by the
# expression-based filter.
# original script: "src/05_10_expression-filter-plots.R"

panel_A <- read_fig_proto("mutation-frequency-hist") +
    theme_fig8() +
    labs(
        y = "number of genes"
    )


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- panel_A

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
