# Build Figure 1

FIGNUM <- 1
VERSION <- 1
FIGFILENAME <- glue("figure_{FIGNUM}_{VERSION}.svg")
FIG_DIMENSIONS <- get_figure_dimensions(2, "medium")


library(patchwork)


#### ---- Figure theme ---- ####

theme_fig1 <- function() {
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


#### ---- A. KRAS mutation lollipop ---- ####

# Panel A.
# The distribution of mutations along the KRAS amino acid sequence.
# original script: "src/90_05_kras-allele-distribution.R"

panel_A <- read_fig_proto("lollipop-kras_2", FIGNUM) &
    theme_fig1() %+replace%
    theme(
        legend.position = c(0.85, 0.85),
        legend.spacing.x = unit(1, "mm"),
        legend.background = element_rect(fill = scales::alpha("white", 0.4),
                                         color = NA),
        legend.margin = margin(1, 1, 1, 1, "mm")
    )

panel_A <- panel_A + labs(tag = "a")


#### ---- B. KRAS allele frequency ---- ####

# Panel B.
# Barplots of the KRAS allele frequency across the cancers.
# original script: "src/90_05_kras-allele-distribution.R"

panel_B <- wrap_plots(
        read_fig_proto("allele_dist_barplot_stackplot", FIGNUM),
        ncol = 2
    ) &
    theme_fig1() %+replace%
    theme(
        axis.title = element_blank(),
        axis.text.y = element_text(hjust = 1)
    )

for (i in 1:4) {
    panel_B[[i]][[2]] <- panel_B[[i]][[2]] *
        theme_fig1() %+replace%
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank()
        )
}

panel_B[[1]][[1]] <- panel_B[[1]][[1]] + labs(tag = "b")


#### ---- C. Mutational signatures probability of causing KRAS allele ---- ####

# Panel C.
# The probability that each allele was caused by each detectable mutational
# signature.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

panel_C <- read_fig_proto("probability-mutsig-caused-allele", FIGNUM) +
    plot_layout(tag_level = "new") &
    theme_fig1() %+replace%
    theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(0, "mm")
    )

panel_C[[1]][[1]] <- panel_C[[1]][[1]] + labs(tag = "c")


#### ---- D. Probability of select signatures causing KRAS allele ---- ####

# Panel D.
# The distribution of the probabilities that alleles were cause by a selection
# of the mutational signatures.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

panel_D <- read_fig_proto("contribution-of-select-signatures", FIGNUM) +
    plot_layout(tag_level = "new") &
    theme_fig1() %+replace%
    theme(
        plot.title = element_text(hjust = 0.5, vjust = 0.7, size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 30, hjust = 1, vjust = 1)
    )

panel_D[[1]] <- panel_D[[1]] + labs(tag = "d")


#### ---- E. Predicting the alleles from the mutational signatures ---- ####

# Panel E.
# The observed vs. predicted KRAS allele frequencies.
# original script: "src/50_10_observed-predicted-kras-alleles.R"

panel_E_proto_list <- read_fig_proto("obs_pred_plot_stats", FIGNUM)
panel_E <- wrap_plots(panel_E_proto_list, nrow = 1) +
    guide_area() +
    plot_layout(guides = "collect",
                widths = c(1, 1, 1, 1, 0.3),
                tag_level = "new") &
    theme_fig1()

panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")


#### ---- F. Predicting the alleles from the mutational signatures ---- ####

# Panel F.
# The observed vs. predicted KRAS allele frequencies.
# original script: "src/50_10_observed-predicted-kras-alleles.R"

panel_F_proto_list <- read_fig_proto("obs_pred_plot_g12_stats", FIGNUM)
panel_F <- wrap_plots(panel_F_proto_list, nrow = 1) +
    guide_area() +
    plot_layout(guides = "collect",
                widths = c(1, 1, 1, 1, 0.3),
                tag_level = "new") &
    theme_fig1()

panel_F[[1]] <- panel_F[[1]] + labs(tag = "f")


#### ---- Figure assembly ---- ####

{
    # ROW 1
    row_1 <- panel_A + panel_B + plot_layout(ncol = 2, widths = c(1, 2, 2))

    # ROW 2
    row_2 <- (panel_C - panel_D) + plot_layout(widths = c(2, 1.2))

    # COMPLETE FIGURE
    full_figure <- row_1 / row_2 / panel_E / panel_F +
        plot_layout(heights = c(1, 1, 0.7, 0.7)) +
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
