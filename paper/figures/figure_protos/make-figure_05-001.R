# Figure 5. Final figure looking at the integration of the comutation and
# genetic dependency results.

FIGNUM <- 5
SUPPLEMENTAL <- FALSE
VERSION <- 1

FIG_DIMENSIONS <- get_figure_dimensions(2, "medium")


theme_fig5 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
    theme_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}

theme_minimal_fig5 <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
    theme_minimal_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#' Special theme for graphs from 'ggraph'.
theme_graph_fig5 <- function() {
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



#### ---- A. Table of comutation and dep. overlaps ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/40_12_overlap-synlet-comutation-hits.R"


panel_A <- read_fig_proto("comut_dep_overlap_tbl.rds", FIGNUM) %>%
    mutate(
        genetic_interaction = str_remove(genetic_interaction, "\\ncomut\\."),
        adj_p_value = scales::scientific(adj_p_value, digits = 3),
        p_val = scales::scientific(p_val, digits = 3)
    )

colnames(panel_A) <- c(
    "Cancer", "Allele", "Gene",
    "Dependency\ncomparison", "Adj. p-value",
    "Comutation\ninteraction", "P-value"
)

core_fontface <- matrix("plain", nrow(panel_A), ncol(panel_A))
core_fontface[, 1:2] <- "bold"
core_fontface[, 3] <- "italic"

ttheme_panel_A <- gridExtra::ttheme_default(
    colhead = list(
        fg_params = list(
            fontsize = 5,
            fontface = "plain",
            fontfamily = "Arial"
        ),
        bg_params = list(fill = "white"),
        padding = unit(c(1, 1), "mm")
    ),
    core = list(
        fg_params = list(
            fontsize = 5,
            fontface = core_fontface,
            fontfamily = "Arial"
        ),
        bg_params = list(fill = "white"),
        padding = unit(c(1, 3), "mm")
    )
)

# Make into a table grob
panel_A <- gridExtra::tableGrob(panel_A, rows = NULL, theme = ttheme_panel_A)

# Add lines
panel_A <- gtable::gtable_add_grob(
    panel_A,
    grobs = grid::segmentsGrob(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0, "npc"), y1 = unit(0, "npc")
    ),
    t = 1, b = 1, l = 1, r = ncol(panel_A)
)
panel_A <- gtable::gtable_add_grob(
    panel_A,
    grobs = grid::segmentsGrob(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0, "npc"), y1 = unit(0, "npc")
    ),
    t = nrow(panel_A), b = nrow(panel_A), l = 1, r = ncol(panel_A)
)

# Wrap for patchwork
panel_A <- wrap_elements(plot = panel_A) +
    labs(tag = "a") +
    theme_fig5()


#### ---- B. Lollipop of STK11 for G12C vs rest ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/40_20_comut-dependency-genes-ppi-connectivity.R"

panel_B <- read_fig_proto("stk11_lollipop_patch", FIGNUM)
panel_B[[1]] <- panel_B[[1]] +
    theme_fig5() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, -2, 0, "mm"),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.position = "bottom"
    ) +
    labs(tag = "b")
panel_B[[2]] <- panel_B[[2]] +
    theme_minimal_fig5() +
    theme(
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.position = "bottom",
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(-2, 0, -2, 0, "mm")
    )
panel_B[[3]] <- panel_B[[3]] +
    theme_fig5() +
    theme(
        plot.margin = margin(-2, 0, 0, 0, "mm"),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.position = "bottom"
    )


#### ---- C. Connectivity of hits ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/40_20_comut-dependency-genes-ppi-connectivity.R"

panel_C <- read_fig_proto("comut_dep_connectivity_bars.rds", FIGNUM) +
    theme_fig5() +
    theme(
        legend.title = element_text(size = 5, hjust = 0),
        legend.key.size = unit(4, "mm"),
        axis.title.x = element_blank()
    ) +
    labs(
        tag = "c"
    )


#### ---- D. Comparing PPI subnets in COAD ---- ####
# A subnetwork of the PPIN for COAD showing the comut. and dep. results for
# all alleles.
# original script: "src/40_15_comparing-COAD-allele-subnetworks.R"

panel_D <- read_fig_proto("coad_overlap_comparison_plot.rds", FIGNUM) +
    theme_graph_fig5() +
    theme(
        legend.position = c(0.8, 0.92),
        legend.title = element_blank()
    ) +
    labs(
        color = "allele",
        tag = "d"
    )



#### ---- Figure assembly ---- ####

{
    set.seed(0)
    # COMPLETE FIGURE
    full_figure <- (
            ((panel_A | panel_B) + plot_layout(widths = c(2, 3))) /
            ((panel_C | panel_D) + plot_layout(widths = c(1, 3)))
        ) +
        plot_layout(heights = c(1, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
