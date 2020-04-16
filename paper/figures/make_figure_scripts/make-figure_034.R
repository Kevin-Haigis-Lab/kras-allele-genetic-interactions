# Figure 034. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 34

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(1, "tall")
FIG_DIMENSIONS$width <- 120

theme_fig34 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. Genetic dependency heatmap ---- ####
# A heatmap of the genes found to have allele-specific genetic dependencies.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

pre_panel_A <- read_fig_proto("PAAD_CRISPR_manhattan_ward.D2_pheatmap")[[4]]

pre_panel_A_main <- gtable::gtable_filter(pre_panel_A,
                                          "legend",
                                          invert = TRUE)

panel_A <- wrap_elements(plot = pre_panel_A_main) *
    theme_fig34(tag_margin = margin(0, -10, 0, 6, "mm")) +
    theme(
        plot.margin = margin(-20, -8, -1, -6, "mm")
    ) +
    labs(tag = "a")



prep_pheatmap_legend <- function(name) {
    read_fig_proto(name) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_fig34() +
        theme(
            plot.margin = margin(0, 3, 0, 3, "mm"),
            plot.title = element_text(size = 6, family = "Arial"),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 6, family = "Arial"),
            plot.background = element_rect(fill = NA, color = NA),
            panel.background = element_rect(fill = NA, color = NA)
        )
}

prep_pheatmap_colorbar <- function(name) {
    read_fig_proto(name) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_fig34() +
            theme(
                plot.margin = margin(0, 3, 0, 3, "mm"),
                plot.title = element_text(size = 6, family = "Arial"),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 6, family = "Arial"),
                plot.background = element_rect(fill = NA, color = NA),
                panel.background = element_rect(fill = NA, color = NA)
            )
}

panel_A_legend1 <- prep_pheatmap_colorbar(
        "PAAD_CRISPR_manhattan_ward.D2_pheatmap_heatpal"
    ) +
    labs(title = "scaled\ndep. score")

panel_A_legend2 <- prep_pheatmap_legend(
        "PAAD_CRISPR_manhattan_ward.D2_pheatmap_allelepal"
    ) +
    labs(title = "allele")

panel_A_legend3 <- prep_pheatmap_legend(
        "PAAD_CRISPR_manhattan_ward.D2_pheatmap_clusterpal"
    ) +
    labs(title = "cluster")

panel_A_legend <- (
        panel_A_legend1 | panel_A_legend2 | panel_A_legend3
    ) +
        plot_layout(widths = c(1, 1, 1))

panel_A_legend <- wrap_elements(full = panel_A_legend) +
    theme_fig34()


#### ---- B. Genetic dependency box-plots ---- ####
# Box-plots of example genes found to have allele-specific genetic dependencies.
# original script: "src/10_11_syn-let_heatmaps-boxplots.R"

panel_B_files <- c(
    "PAAD-CEP350",
    "PAAD-EGLN2",
    "PAAD-JUN",
    "PAAD-NUMB",
    "PAAD-TMED2"
)

panel_B_plots <- as.list(rep(NA, length(panel_B_files)))
names(panel_B_plots) <- panel_B_files
for (f in panel_B_files) {

    panel_B_plots[[f]] <- read_fig_proto(f) +
        theme_fig34(tag_margin = margin(0, 0, 0, -2, "mm")) %+replace%
        theme(
            plot.title = element_text(size = 6, face = "bold"),
            axis.title.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(0, 0, 3, 0, "mm")
        ) +
        labs(
            title = str_remove(file_sans_ext(f), "PAAD-")
        )

}

panel_B <- wrap_plots(panel_B_plots, ncol = 1)
panel_B[[1]] <- panel_B[[1]] + labs(tag = "b")



#### ---- Figure assembly ---- ####


col1_layout <- "
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    ##BBB##
"

{
    set.seed(0)


    column1 <- (panel_A + panel_A_legend) + plot_layout(design = col1_layout)

    # COMPLETE FIGURE
    full_figure <- (
        column1 |
        (
            (panel_B / plot_spacer()) +
            plot_layout(heights = c(5, 5, 5, 5, 5, 2))
        )
    ) +
        plot_layout(widths = c(7, 2))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
