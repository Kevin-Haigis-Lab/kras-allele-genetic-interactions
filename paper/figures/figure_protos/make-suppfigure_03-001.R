
# Supplemental Figure 2. Distribution of mutational signatures across alleles

FIGNUM <- 3
SUPPLEMENTAL <- TRUE
VERSION <- 1
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")

FIG_DIMENSIONS$height <- 100  # specific change for this figure


theme_figS3 <- function() {
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

panel_A <- read_fig_proto("mutation-frequency-hist",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_figS3()



#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- panel_A +
        plot_annotation(
            title = glue("Supp Figure {FIGNUM}"),
            theme = theme(
                plot.title = element_text(size = 10,
                                          family = "Arial",
                                          hjust = 0)
            )
        )

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
