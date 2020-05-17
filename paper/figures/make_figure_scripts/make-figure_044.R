# Figure 044. Observed vs. predicted KRAS allele frequencies with all alleles.

FIGNUM <- 44

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 165


theme_fig44 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
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
                        "_predict-ALL-allele-freq_scatter")

codon_pal <- codon_palette[names(codon_palette) != "Other"]


get_max_val <- function(gg) {
    data <- ggplot_build(read_fig_proto(gg))$plot$data
    return(max(c(data$upper_ci, data$observed_allele_frequency)))
}

panel_A_proto_list <- map(panel_A_plots, function(x) {

    max_val <- get_max_val(x)
    min_val <- (-0.1 / 0.45) * max_val

    p <- read_fig_proto(x) +
        scale_x_continuous(limits = c(min_val, max_val),
                           expand = expansion(mult = c(0, 0.03))) +
        scale_y_continuous(limits = c(min_val, max_val),
                           expand = expansion(mult = c(0, 0.03)),
                           labels = function(x) {ifelse(x == 0, "", x)}) +
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

panel_A <- wrap_plots(panel_A_proto_list, nrow = 2, guides = "collect") &
    theme_fig44() %+replace%
    theme(
        plot.title = element_text(size = 7, vjust = 2, face = "bold"),
        plot.subtitle = element_markdown(hjust = 0, vjust = 0.3),
        legend.position = "right",
        legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(1, "mm"),
        legend.key.width = unit(1, "mm"),
        legend.spacing.y = unit(3, "mm"),
        plot.margin = margin(3, 3, 3, 3, "mm")
    )

for (i in seq(2, 4)) {
    panel_A[[i]] <- panel_A[[i]] + labs(y = NULL)
}

panel_A <- (panel_A | guide_area()) +
    plot_layout(widths = c(20, 1))


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- panel_A +
        plot_layout()

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
