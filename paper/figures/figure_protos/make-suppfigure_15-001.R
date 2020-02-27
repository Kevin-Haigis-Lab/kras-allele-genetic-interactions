# Supplemental Figure 15. "Increased expression of WDR26 has previously been
# implicated in driving breast cancer by serving as a scaffolding protein in
# the PI3K-Akt pathway, though there does not appear to be a link between RNA
# levels and genetic dependency in the COAD cell lines"

FIGNUM <- 15
SUPPLEMENTAL <- TRUE
VERSION <- 1

#> SET THE FIGURE DIMENSIONS
# FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS <- list(width = 110, height = 100)


theme_figS15 <- function() {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


#### ---- A. WDR26 dependency by mRNA expression ---- ####
# A scatter plot of WDR26 dependency by mRNA expression in COAD
# original script: "src/10_57_coad_depmap_wdr26.R"

panel_A <- read_fig_proto("coad_depmap_wdr26-rna-v-dep.svg",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_figS15() +
    labs(tag = "a")


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- panel_A +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
