# Figure 043. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 43

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 50


theme_fig43 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/50_12_observed-predicted-kras-alleles_v3.R"

panel_A_plots <- paste0(c("COAD", "LUAD", "MM", "PAAD"),
                        "_predict-allele-freq_scatter.svg")

codon_pal <- codon_palette[names(codon_palette) != "Other"]

panel_A_proto_list <- map(panel_A_plots, function(x) {
    p <- read_fig_proto(x) +
        scale_color_manual(
            values = codon_pal,
            drop = FALSE,
            guide = guide_legend(
                order = 10,
                title.theme = element_markdown(hjust = 0.5, vjust = 0,
                                               face = "bold", size = 6,
                                               family = "Arial"),
                label.vjust = 0.5,
                keyheight = unit(4, "mm"),
                label.position = "right"
            )
        ) +
        scale_shape_manual(
            values = c(17, 16),
            drop = FALSE,
            guide = guide_legend(
                order = 20,
                title.position = "top",
                title.theme = element_markdown(hjust = 0.5, vjust = 0,
                                               face = "bold", size = 6,
                                               family = "Arial"),
                label.position = "right",
                label.vjust = 0.5,
                keyheight = unit(3, "mm")
            )
        )
    return(p)
})

panel_A <- wrap_plots(panel_A_proto_list, nrow = 1, guides = "collect") &
    theme_fig43() %+replace%
    theme(
        plot.title = element_text(size = 7, vjust = 2, face = "bold"),
        plot.subtitle = element_markdown(hjust = 0, vjust = 0.3),
        legend.position = "right",
        legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(1, "mm"),
        legend.key.width = unit(1, "mm"),
        legend.spacing.y = unit(3, "mm"),
        plot.margin = margin(0, 1, 0, 1, "mm")
    )

for (i in seq(2, 4)) {
    panel_A[[i]] <- panel_A[[i]] + labs(y = NULL)
}

panel_A <- (panel_A | guide_area()) +
    plot_layout(widths = c(4, 4, 4, 4, 1))


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
