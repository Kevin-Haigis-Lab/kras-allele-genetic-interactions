
# Supplemental Figure 2. Distribution of mutational signatures across alleles

FIGNUM <- 2
SUPPLEMENTAL <- TRUE
VERSION <- 1
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")




theme_figS2 <- function() {
    theme_bw(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5, hjust = 1),
        axis.text.x = element_text(size = 5, vjust = 1),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm")),
        strip.background = element_blank()
    )
}


#### ---- A. Distirubiton of mutational signatures in each sample ---- ####

# The distribution of mutational signature levels in each sample (barplot).
# original script: "src/50_20_mutsignatures-distributions.R"

panel_A <- read_fig_proto("signature-level-per-sample",
                          FIGNUM, supp = SUPPLEMENTAL) +
    scale_fill_manual(
        values = mutsig_descrpt_pal,
        guide = guide_legend(
            nrow = 1,
            label.position = "top",
            label.hjust = 0.5,
            label.vjust = -9
        )
    ) +
    labs(
        y = "signature level",
        tag = "a"
    ) &
    theme_figS2() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(4, "mm")
    )



#### ---- B. Distirubiton of mutational signature levels ---- ####

# The distribution of mutational signature levels in each sample (boxplot).
# original script: "src/50_20_mutsignatures-distributions.R"

panel_B <- read_fig_proto("signature-level-boxplots_with0",
                          FIGNUM, supp = SUPPLEMENTAL) +
    labs(
        tag = "b",
        y = "signature level"
    ) &
    theme_figS2() +
    theme(
        legend.position = "none"
    )

#### ---- C. Distirubiton of mutational signatures by allele ---- ####

# The distribution of mutational signature levels in samples with each allele.
# original script: "src/50_20_mutsignatures-distributions.R"

panel_C <- read_fig_proto("mutational-signatures-distribution-by-allele",
                          FIGNUM, supp = SUPPLEMENTAL)

for (i in seq(1, length(panel_C))) {
    panel_C[[i]] <- panel_C[[i]] +
        labs(
            y = "level"
        )
}

panel_C[[1]] <- panel_C[[1]] + labs(tag = "c")

panel_C <- panel_C %>%
    wrap_plots(nrow = 2) &
    theme_figS2() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank()
    )


#### ---- D. Levels of clock vs. non-clock mutational signatures ---- ####

# The levels of clock and non-clock mutational signatures per cancer.
# original script: "src/50_20_mutsignatures-distributions.R"

panel_D <- read_fig_proto("clock-signature-plots",
                          FIGNUM, supp = SUPPLEMENTAL)

for (i in seq(1, length(panel_D))) {
    panel_D[[i]] <- panel_D[[i]] +
        labs(
            y = "level"
        )
}

panel_D[[1]] <- panel_D[[1]] + labs(tag = "d")

panel_D <- panel_D %>%
    wrap_plots(nrow = 2) &
    theme_figS2() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank()
    )


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (panel_A | panel_B) / (panel_C | panel_D)  / guide_area() +
        plot_layout(heights = c(15, 10, 1), guides = "collect") +
        plot_annotation(
            title = glue("Supp Figure {FIGNUM}"),
            theme = theme(
                plot.title = element_text(size = 10,
                                          family = "Arial",
                                          hjust = 0)
            )
        )

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS,
        file_fmt = "jpeg"
    )
    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
