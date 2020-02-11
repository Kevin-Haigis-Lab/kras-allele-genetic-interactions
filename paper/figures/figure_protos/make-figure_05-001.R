# Figure 5. Final figure looking at the integration of the comutation and
# genetic dependency results.

library(gridExtra)
library(gtable)
library(grid)

FIGNUM <- 5
SUPPLEMENTAL <- FALSE
VERSION <- 1

FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig5 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


#### ---- A. Table of comutation and dep. overlaps ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/40_12_overlap-synlet-comutation-hits.R"


panel_A <- read_fig_proto("comut_dep_overlap_tbl", FIGNUM) %>%
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

ttheme_panel_A <- gridExtra::ttheme_default(
    core = list(
        fg_params = list(fontsize = 5)),
    colhead = list(
        fg_params = list(fontsize = 5, fontface = "bold"),
        bg_params = list(fill = "grey80")
    )
)

panel_A <- tableGrob(panel_A, rows = NULL, theme = ttheme_panel_A)

panel_A <- gtable::gtable_add_grob(
    panel_A,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 2, b = nrow(panel_A), l = 1, r = ncol(panel_A)
)

panel_A <- wrap_elements(panel = panel_A)


#### ---- B. Connectivity of hits ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/40_12_overlap-synlet-comutation-hits.R"

panel_B <- plot_spacer()


#### ---- C. Comparing PPI subnets in COAD ---- ####
# A table of the genes found to comutate and show differential dependency
# with an allele.
# original script: "src/40_12_overlap-synlet-comutation-hits.R"

panel_C <- plot_spacer()



#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (
            ((panel_A | panel_B) + plot_layout(widths = c(1, 1))) /
            panel_C
        ) +
        plot_layout(heights = c(1, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
