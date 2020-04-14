# __FIGURE_NAME__. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- __FIGURE_NUM__

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig__THEME_SUFFIX__ <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. TITLE ---- ####
# A BRIEF DESCRIPTION OF THE PANEL.
# original script: "src/##_##_ORIGINAL-SCRIPT.R"

panel_A <- #> BEGIN PANEL A HERE


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- panel_A +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
