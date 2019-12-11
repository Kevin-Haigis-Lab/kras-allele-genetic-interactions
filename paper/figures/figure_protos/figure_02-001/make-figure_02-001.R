# Build Figure 2

FIGNUM <- 2
VERSION <- 1
FIGFILENAME <- glue("figure_{FIGNUM}_{VERSION}.svg")
FIG_DIMENSIONS <- get_figure_dimensions(2, "medium")


library(patchwork)


#### ---- Figure theme ---- ####

theme_fig2 <- function() {
    theme_bw(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5, hjust = 1),
        axis.text.x = element_text(size = 5, vjust = 1),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-3, -3, -3, -3, "mm"))
    )
}

theme_graph_fig2 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        legend.position = "left",
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-3, -3, -3, -3, "mm")),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}


#### ---- A. High-level comutation network for COAD ---- ####

# Panel A.
# The high-level network plot for the comutation graph for COAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_COAD.svg", 2) *
    theme_graph_fig2() %+replace%
    theme(
        legend.spacing.x = unit(1.5, "mm")
    )

panel_A <- guide_area() + (panel_A + labs(tag = "a")) +
    plot_layout(guides = "collect", widths = c(1, 10))


#### ---- B. A priori genes of interest comutation network for COAD ---- ####

# Panel B.
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"

panel_B <- read_fig_proto("goi_overlap_genetic_interactions_network_COAD_allLists", 2) *
    theme_graph_fig2() +
    theme(
        legend.position = "bottom"
    )

panel_B <- panel_B + labs(tag = "b")


#### ---- C. Rainfall plot ---- ####

# Panel C.
# A rainfall plot.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"

library(grImport2)
library(grid)

plot_file_name <- "COAD_G12D_exclusivity_oncostrip_select.svg"
graphs_dir <- file.path("graphs", "20_50_rainfall-plots-select")
file_path <- file.path(graphs_dir, plot_file_name)

oncoplot_pic <- readPicture(file_path)
oncoplot_grob <- gTree(children = gList(pictureGrob(oncoplot_pic)))

panel_C <- ggplot(tibble()) +
    annotation_custom(oncoplot_grob) +
    theme(
        panel.background = element_blank()
    )


#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    # ROW 1
    row_1 <- (panel_A | panel_B) / (panel_C + plot_spacer()) +
        plot_layout(heights = c(2, 1))

    # COMPLETE FIGURE
    full_figure <- row_1 +
        plot_annotation(
            title = glue("Figure {FIGNUM}"),
            theme = theme(
                plot.title = element_text(size = 10, family = "Arial")
            )
        )

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
