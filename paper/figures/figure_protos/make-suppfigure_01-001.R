
# Supplemental Figure 1. Barplots of frequency of all alleles

FIGNUM <- 1
VERSION <- 1

FIG_DIMENSIONS <- get_figure_dimensions(2, "short")

FIG_DIMENSIONS$height <- 100  # specific change for this figure


#### ---- Figure theme ---- ####

theme_figS1 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-3, -3, -3, -3, "mm"))
    )
}



#### ---- A. Distribution of all alleles in each cancer ---- ####

# The distribution of the alleles (no "Other" group) per cancer.
# original script: "src/90_05_kras-allele-distribution.R"

panel_A <- read_fig_proto("allele_dist_barplot_stackplot_all",
                         FIGNUM,
                         supp = TRUE) %>%
    wrap_plots() &
    theme_figS1() %+replace%
    theme(
        axis.title.y = element_blank(),
    )

for (i in 1:4) {
    panel_A[[i]][[2]] <- panel_A[[i]][[2]] *
        theme_figS1() %+replace%
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank()
        )
}

#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- panel_A

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = TRUE,
        dim = FIG_DIMENSIONS
    )
}
