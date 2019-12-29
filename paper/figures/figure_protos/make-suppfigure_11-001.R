# Supplemental Figure 11. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 11
SUPPLEMENTAL <- TRUE
VERSION <- 1

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_figS11 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}

#' Special theme for graphs from 'ggraph'.
theme_graph_figS11 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm")),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "bottom",
        legend.title = element_text(size = 5, hjust = 0.5),
        legend.text = element_text(size = 5, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}



#### ---- A. PPI of Myc ---- ####
# The PPI of Myc (TF) - an enriched function in the comutation network of
# G12C in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_A <- read_fig_proto("LUAD_G12C_Transcription-Factor-PPIs_MYC",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_graph_figS11() +
    labs(tag = "a")


#### ---- B. PPI of focal adhesions ---- ####
# The PPI of focal adhesions - an enriched function in the comutation network of
# G12D in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_B <- read_fig_proto("LUAD_G12D_KEGG-2019-Human_Focal-adhesion",
                          FIGNUM, supp = SUPPLEMENTAL) +
    theme_graph_figS11() +
    labs(tag = "b")


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- panel_A / panel_B +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
