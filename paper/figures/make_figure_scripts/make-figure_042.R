# Figure 042. #> Mutational signature main figure.
# (Based off of Fig 28.)

FIGNUM <- 42

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 160

theme_fig42 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}



#### ---- A. Distribution of KRAS alleles ---- ####
# The distribution of KRAS alleles across cancers and codon.
# original script: "src/90_05_kras-allele-distribution.R"

panel_A <- read_fig_proto("allele_dist_dotplot") +
    scale_size_continuous(labels = scales::label_percent(accuracy = 1)) +
    theme_fig42(tag_margin = margin(-1, -1, -1, -5.5, "mm")) +
    theme(
        axis.ticks = element_blank(),
        axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.position = "left",
        legend.text = element_text(size = 6),
        legend.title = element_markdown(size = 6),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold")
    ) +
    labs(tag = "a",
         size = "percent of<br>*KRAS* mutations")


#### ---- B. Percent of samples with KRAS mutation ---- ####
# Percent of samples with a KRAS mutation per cancer.
# original script: "src/90_05_kras-allele-distribution.R"

panel_B <- read_fig_proto("cancer_freq_kras_mut_column") +
    scale_x_continuous(
        expand = expansion(add = c(0, 0.02)),
        breaks = c(0.2, 0.4, 0.6, 0.8),
        labels = function(x) { paste0(round(x * 100), "%") }
    ) +
    theme_fig42() +
    theme(
        axis.title.x = element_markdown(),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank()
    ) +
    labs(tag = "b")



#### ---- C. Distirubiton of mutational signatures by allele ---- ####
# The distribution of mutational signature levels in samples with each allele.
# original script: "src/50_20_mutsignatures-distributions.R"

mutsig_barplot_widths <- c(7, 4, 5, 5)

style_mutsig_prob_barplots <- function(plt, i, tag = NULL, y = NULL) {
    themed_plt <- plt +
        scale_fill_manual(
            values = mutsig_descrpt_pal,
            guide = guide_legend(
                nrow = 1,
                title.position = "left",
                title.vjust = 0.2,
                label.position = "top",
                label.hjust = 0.5,
                label.vjust = -4.5
            )
        ) +
        theme_fig42(tag_margin = margin(-1, -1, -1, -9, "mm")) +
        theme(
            plot.title = element_text(vjust = 0.5, size = 7, face = "bold"),
            plot.margin = margin(1, 1, 1, 1, "mm"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 6),
            legend.position = "none",
            legend.title = element_text(size = 6, face = "bold"),
            legend.text = element_text(size = 6),
            legend.key.size = unit(3, "mm"),
            legend.spacing.x = unit(1, "mm"),
            legend.spacing.y = unit(0, "mm")
        )
    if (i == 1) {
        themed_plt <- themed_plt +
            theme(axis.text.y = element_text(size = 6)) +
            labs(tag = tag, y = y)
    } else {
        themed_plt <- themed_plt +
            theme(axis.text.y = element_blank()) +
            labs(y = " ")
    }
    return(themed_plt)
}

panel_C_plots <- paste0(c("COAD", "LUAD", "MM", "PAAD"),
                        "_mutational-signatures-distribution-by-allele")

panel_C <- lapply(panel_C_plots, read_fig_proto) %>%
    imap(style_mutsig_prob_barplots, tag = "c", y = "signature composition") %>%
    wrap_plots(nrow = 1, widths = mutsig_barplot_widths)


#### ---- D. Mutational signatures probability of causing KRAS allele ---- ####
# The probability that each allele was caused by each detectable mutational
# signature.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

panel_D <- read_fig_proto("probability-mutsig-caused-allele_barplot-list") %>%
    imap(style_mutsig_prob_barplots, tag = "d", y = "probability of causing allele")


pull_signatures_from_panel_D <- function(x) {
    unique(ggplot_build(x)$plot$data$description)
}

signatures <- map(panel_D, pull_signatures_from_panel_D) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    as.character()

panel_D_legend <- custom_label_legend(
        signatures,
        y_value = "signature",
        family = "Arial",
        size = 1.8,
        label.padding = unit(1, "mm"),
        label.size = unit(0, "mm"),
        hjust = 0.5
    ) +
    scale_fill_manual(values = mutsig_descrpt_pal) +
    theme(
        legend.position = "none",
        plot.margin = margin(-10, 0, -10, 0, "mm"),
        axis.text.y = element_text(size = 6, face = "bold")
    )

panel_D <- wrap_plots(panel_D, nrow = 1, widths = mutsig_barplot_widths)


#### ---- E. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/50_12_observed-predicted-kras-alleles_v3.R"

panel_E_plots <- paste0(c("COAD", "LUAD", "MM", "PAAD"),
                        "_predict-allele-freq_scatter.svg")

panel_E_proto_list <- lapply(panel_E_plots, read_fig_proto)
panel_E <- wrap_plots(panel_E_proto_list, nrow = 1, guides = "collect") &
    theme_fig42() %+replace%
    theme(
        plot.title = element_text(size = 7, face = "bold"),
        legend.title = element_markdown(vjust = 0.5, hjust = 0.5,
                                        face = "bold", size = 6,
                                        family = "Arial"),
        legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(1, "mm"),
        legend.key.width = unit(1, "mm"),
        legend.spacing.y = unit(3, "mm"),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )

for (i in seq(2, 4)) {
    panel_E[[i]] <- panel_E[[i]] + labs(y = "")
}

for (i in seq(2, 4)) {
    panel_E[[i]] <- panel_E[[i]] + theme(legend.position = "none")
}

panel_E[[1]] <- panel_E[[1]] + labs(tag = "e")

panel_E <- (panel_E | guide_area()) +
    plot_layout(widths = c(10, 10, 10, 10, 3))


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    row_1 <- (panel_A | panel_B) +
        plot_layout(widths = c(10, 3))

    panel_D_legend_sp <- (plot_spacer() | panel_D_legend | plot_spacer()) +
        plot_layout(widths = c(1, 6, 1))

    row_2 <- (panel_C / panel_D / panel_D_legend_sp) +
        plot_layout(heights = c(6, 6, 1))

    row_2 <- wrap_elements(full = row_2)

    full_figure <- (row_1 / row_2 / wrap_elements(full = panel_E)) +
        plot_layout(heights = c(2, 5, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
