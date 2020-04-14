# Figure 027. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 27

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")
FIG_DIMENSIONS$height <- 120

theme_fig27 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


theme_graph_fig16 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "bottom",
        legend.title = element_text(size = 5, hjust = 0.5),
        legend.text = element_text(size = 5, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}



#### ---- A. PPI of focal adhesions ---- ####
# The PPI of focal adhesions - an enriched function in the comutation network of
# G12D in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_A <- read_fig_proto("LUAD_G12C_Transcription-Factor-PPIs_MYC") +
    theme_graph_fig16() +
    scale_color_manual(
        drop = FALSE,
        values = c(comut_updown_pal,
                   "none" = "grey70",
                   "in_geneset" = "grey40"),
        label = function(x) { str_replace_all(x, "_", " ") },
        guide = guide_legend(
            title = "comutation interaction",
            nrow = 2,
            title.theme = element_text(size = 6, face = "bold", hjust = 0.5),
            title.position = "top",
            label.theme = element_text(size = 6, hjust = 0),
            label.position = "right",
            label.hjust = 0
        )
    ) +
    theme(
        legend.position = c(0.25, 0.95),
        legend.key.height = unit(3, "mm")
    ) +
    labs(tag = "a")



#### ---- B. PPI of Myc ---- ####
# The PPI of Myc (TF) - an enriched function in the comutation network of
# G12C in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_B <- read_fig_proto("LUAD_G12D_KEGG-2019-Human_Focal-adhesion") +
    theme_graph_fig16() +
    theme(
        legend.position = "none"
    ) +
    labs(tag = "b")


#### ---- Figure assembly ---- ####

{
    set.seed(2)

    figure_27_design <- "
        111111111222
        111111111222
        111111111###
        111111111###
    "

    # COMPLETE FIGURE
    full_figure <- (panel_A  + panel_B) +
        plot_layout(design = figure_27_design)

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
