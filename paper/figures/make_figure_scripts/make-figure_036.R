# Figure 036. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 36

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig36 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}

# Special theme for graphs from 'ggraph'.
theme_graph_fig36 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm")),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "right",
        legend.title = element_text(size = 5, hjust = 0.5),
        legend.text = element_text(size = 5, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}


#### ---- A. PPIN for LUAD G12C ---- ####
# Annotated PPIN for LUAD G12C.
# original script: "src/40_16_comparing-LUAD-allele-subnetworks.R"

panel_A <- read_fig_proto("luad-G12C_overlap_comparison_plot") +
    theme_graph_fig36() +
    theme(
        legend.position = "none"
    )


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
