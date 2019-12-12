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


adjust_oncoplot_theme <- function(pw) {
    # Top bar plot
    pw[[1]] <- pw[[1]] +
        theme(
            axis.text.y = element_text(size = 5, hjust = 1),
            plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm"))
        )

    # Middle main tile plot
    pw[[2]] <- pw[[2]] +
        theme(axis.text.y = element_text(size = 6, hjust = 1))

    # Top bar plot
    pw[[3]] <- pw[[3]] +
        theme(axis.text.x = element_text(size = 5, vjust = 1))

    return(pw)
}


#### ---- A. High-level comutation network for COAD ---- ####

# Panel A.
# The high-level network plot for the comutation graph for COAD.
# original script: "src/20_40_highlivel-genetic-interactions.R"

panel_A <- read_fig_proto("genetic_interaction_network_COAD", 2) +
    theme_graph_fig2() %+replace%
    theme(
        legend.spacing.x = unit(1, "mm")
    ) +
    labs(tag = "a")


#### ---- B. A priori genes of interest comutation network for COAD ---- ####

# Panel B.
# The subset from the high-level network for genes known to be related to KRAS.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"

panel_B <- read_fig_proto("goi_overlap_genetic_interactions_network_COAD_allLists", 2) *
    theme_graph_fig2() +
    theme(
        legend.position = "bottom"
    ) +
    labs(tag = "b")


#### ---- C. Lollipop ---- ####

# Panel C.
# A lollipop plot.
# original script: "src/20_43_apriori-lists-genetic-interactions.R"



#### ---- D. Oncoplot ---- ####

# Panel D.
# A rainfall plot.
# original script: "src/20_50_rainfall-plots.R"

panel_D <- read_fig_proto("COAD_G12D_comutation_oncostrip_select", FIGNUM)
panel_D <- adjust_oncoplot_theme(panel_D)
panel_D[[1]] <- panel_D[[1]] + labs(tag = "d") +
    theme()


#### ---- E. Oncoplot ---- ####

# Panel E.
# A rainfall plot.
# original script: "src/20_50_rainfall-plots.R"

panel_E <- read_fig_proto("COAD_G12D_exclusivity_oncostrip_select", FIGNUM)
panel_E <- adjust_oncoplot_theme(panel_E)
panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")


#### ---- Figure assembly ---- ####

{
    set.seed(0)  # Because the graph-plotting algorithm is stochastic.

    row_1 <- (panel_A | panel_B | plot_spacer()) +
        plot_layout(widths = 1)

    row_2 <- (plot_spacer() | panel_D | panel_E) +
              plot_layout(widths = c(1, 1000, 1000))

    row_3 <- plot_spacer()

    # COMPLETE FIGURE
    full_figure <- (row_1) / (row_2)  / (row_3) +
        plot_layout(heights = c(2, 1, 2)) +
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
