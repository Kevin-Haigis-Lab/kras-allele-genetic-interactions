# Figure 028. #> Mutational signature main figure.
# (Based off of sketch from PJP.)

FIGNUM <- 28

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "short")


theme_fig28 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. KRAS hot-spot mutation frequency barplot ---- ####
# The distribution of mutations at the 4 hotspots.
# original script: "src/90_05_kras-allele-distribution.R"

codons_to_label <- c(12, 13, 61, 146)

panel_A <- read_fig_proto("lollipop-kras_hotspot-only") +
    theme_fig28() +
    theme(
        panel.grid.major.x = element_blank(),
        axis.title.y = element_markdown(),
        legend.position = c(0.80, 0.82),
        legend.background = element_blank(),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6),
        plot.margin = margin(0, 2, 0, 0, "mm")
    ) +
    labs(tag = "a",
         y = "number of samples (*log*<sub>10</sub>)")


#### ---- B. Distribution of KRAS alleles ---- ####
# The distribution of KRAS alleles across cancers and codon.
# original script: "src/90_05_kras-allele-distribution.R"

panel_B <- read_fig_proto("allele_dist_dotplot") +
    theme_fig28() +
    theme(
        axis.ticks = element_blank(),
        axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold")
    ) +
    labs(tag = "b")


#### ---- C. Percent of samples with KRAS mutation ---- ####
# Percent of samples with a KRAS mutation per cancer.
# original script: "src/90_05_kras-allele-distribution.R"

panel_C <- read_fig_proto("cancer_freq_kras_mut_column") +
    scale_x_continuous(
        expand = expansion(add = c(0, 0.02)),
        breaks = c(0.2, 0.4, 0.6, 0.8),
        labels = function(x) { paste0(round(x * 100), "%") }
    ) +
    theme_fig28() +
    theme(
        axis.title.x = element_markdown(),
        axis.text.x = element_text(hjust = 0, vjust = 0.5, angle = -90),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank()
    ) +
    labs(tag = "c")


#### ---- D. Mutational signature demonstration ---- ####
# Demonstration of mutational signatures decomposing a mutational spectrum.
# original script: "src/50_60_example-mutational-signature-spectra.R"

mutational_spectrum_theme <- function(...) {
    theme_fig28(...) %+replace%
    theme(
        plot.title = element_text(size = 7, face = "bold", vjust = 1.4),
        axis.text.x = element_text(size = 5, angle = -90,
                                   hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        legend.position = "none",
        panel.spacing.x = unit(0, "mm"),
        panel.grid.major.x = element_blank()
    )
}

panel_D_1 <- read_fig_proto("coad_mut_spectrum") +
    mutational_spectrum_theme() +
    labs(title = "COAD mutational spectrum",
         tag = "d")


mutational_spectrum_theme_minimal <- function(...) {
    mutational_spectrum_theme() %+replace%
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", size = 5),
    )
}


read_panel_D_subpanel <- function(sig) {
    glue("sig{sig}_example_mutsig_spectra") %>%
        as.character() %>%
        read_fig_proto() +
        mutational_spectrum_theme_minimal() +
        labs(title = glue("Sig. {sig}"))
}
panel_D_subpanels <- c("1", "4", "5", "8", "9", "18") %>%
    map(read_panel_D_subpanel) %>%
    wrap_plots(nrow = 2)

panel_D <- (panel_D_1 / panel_D_subpanels)



#### ---- E. Mutational signatures probability of causing KRAS allele ---- ####
# The probability that each allele was caused by each detectable mutational
# signature.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

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
        theme_fig28() +
        theme(
            plot.title = element_text(vjust = 0.5, size = 7, face = "bold"),
            plot.margin = margin(1, 1, 1, 1, "mm"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 5.8),
            legend.position = "none",
            legend.title = element_text(size = 6, face = "bold"),
            legend.text = element_text(size = 6),
            legend.key.size = unit(3, "mm"),
            legend.spacing.x = unit(1, "mm"),
            legend.spacing.y = unit(0, "mm")
        )
    if (i == 1) {
        themed_plt <- themed_plt + labs(tag = tag)
    }
    if (i %% 2 == 1) {
        themed_plt <- themed_plt + labs(y = y)
    } else {
        themed_plt <- themed_plt + labs(y = NULL)
    }
    return(themed_plt)
}


panel_E <- read_fig_proto("probability-mutsig-caused-allele_barplot-list") %>%
    imap(style_mutsig_prob_barplots, tag = "e", y = "probability")

panel_E_design <- "
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    111111122222
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    333333444444
    555555555555
"

pull_signatures_from_panel_E <- function(x) {
    unique(ggplot_build(x)$plot$data$description)
}

signatures <- map(panel_E, pull_signatures_from_panel_E) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    as.character()

panel_E_legend <- custom_label_legend(
        signatures,
        y_value = "signature",
        family = "Arial", size = 1.8,
        label.padding = unit(1, "mm"),
        label.size = unit(0, "mm"),
        hjust = 0.5
    ) +
    scale_fill_manual(values = mutsig_descrpt_pal) +
    theme(
        legend.position = "none",
        plot.margin = margin(-5, 0, 0, 0, "mm")
    )

panel_E <- wrap_plots(panel_E) +
    plot_layout(design = panel_E_design)


#### ---- F. Predicted vs Observed allele frequency ---- ####
# Predicted vs. observed frequency of KRAS alleles in each cancer.
# original script: "src/50_10_observed-predicted-kras-alleles.R"

panel_F_plots <- paste0(c("COAD", "LUAD", "MM", "PAAD"),
                        "_predict-allele-freq_scatter.svg")

panel_F_proto_list <- lapply(panel_F_plots, read_fig_proto)
panel_F <- wrap_plots(panel_F_proto_list, nrow = 1, guides = "collect") &
    theme_fig28() %+replace%
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

for (i in seq(2, 4)) panel_F[[i]] <- panel_F[[i]] + labs(y = "")

panel_F[[1]] <- panel_F[[1]] + labs(tag = "f")

panel_F <- (panel_F | guide_area()) +
    plot_layout(widths = c(10, 10, 10, 10, 3))


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    row_1 <- (panel_A | panel_B | panel_C) +
        plot_layout(widths = c(4, 10, 3))

    panel_E_legend_spacer <- (plot_spacer() | panel_E_legend | plot_spacer()) +
        plot_layout(widths = c(1, 20, 1))

    panel_E_group <- (panel_E / panel_E_legend) +
            plot_layout(heights = c(15, 1))

    row_2 <- (panel_D | wrap_elements(full = panel_E_group)) +
        plot_layout(widths = c(8, 10))

    full_figure <- (row_1 / row_2 / wrap_elements(full = panel_F)) +
        plot_layout(heights = c(2, 3, 2))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
