# Figure 021. PPIN from integration analysis for LUAD and PAAD.

FIGNUM <- 21

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig21 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


# Special theme for graphs from 'ggraph'.
theme_graph_fig21 <- function() {
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
    theme_graph_fig21() +
    theme(
        legend.position = "none"
    ) +
    labs(
        tag = "a"
    )


#### ---- B. PPIN for LUAD G12V ---- ####
# Annotated PPIN for LUAD G12V.
# original script: "src/40_16_comparing-LUAD-allele-subnetworks.R"

panel_B <- read_fig_proto("luad-G12V_overlap_comparison_plot") +
    theme_graph_fig21() +
    theme(
        legend.position = "none"
    ) +
    labs(
        tag = "b"
    )


#### ---- C. PPIN for PAAD ---- ####
# Annotated PPIN for PAAD.
# original script: "src/40_17_comparing-PAAD-allele-subnetworks.R"

panel_C <- read_fig_proto("paad_overlap_comparison_plot") +
    theme_graph_fig21() +
    theme(
        legend.position = c(0.1, 0.80),
        legend.title = element_blank()
    ) +
    labs(
        color = "allele",
        tag = "c"
    )



#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- (panel_A / panel_B / panel_C) +
        plot_layout(heights = c(3, 2, 2))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
