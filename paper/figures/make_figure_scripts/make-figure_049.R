# Figure 049. #> Associations between mutational signature level and KRAS
#   alleles.

FIGNUM <- 49

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig49 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}

theme_msig_heatmap <- function() {
    theme_fig49() %+replace%
        theme(axis.title = element_blank(),
              legend.title = element_text(),
              strip.background = element_blank(),
              strip.text = element_text(size = 7, face = "bold"),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = 6))
}


theme_msig_boxplot <- function() {
    theme_fig49() %+replace%
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                       size = 6),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 7, face = "bold")
        )
}


#### ---- A. COAD ---- ####
# Associations between KRAS alleles and mutational sigs. in COAD
# original script: "src/50_35_mutational-signature-allele-associations.R"

panel_A1 <- read_fig_proto("COAD_median_allele-sig-associate.svg") +
    theme_msig_heatmap() +
    labs(tag = "a")

panel_A2 <- read_fig_proto("COAD_allele-sig-boxplots") +
    theme_msig_boxplot()



#### ---- B. LUAD ---- ####
# Associations between KRAS alleles and mutational sigs. in LUAD
# original script: "src/50_35_mutational-signature-allele-associations.R"

panel_B1 <- read_fig_proto("LUAD_median_allele-sig-associate.svg") +
    theme_msig_heatmap() +
    labs(tag = "b")

panel_B2 <- read_fig_proto("LUAD_allele-sig-boxplots") +
    theme_msig_boxplot()



#### ---- C. MM ---- ####
# Associations between KRAS alleles and mutational sigs. in MM
# original script: "src/50_35_mutational-signature-allele-associations.R"

panel_C1 <- read_fig_proto("MM_median_allele-sig-associate.svg") +
    theme_msig_heatmap() +
    theme(legend.box = "horizontal") +
    labs(tag = "c")

panel_C2 <- read_fig_proto("MM_allele-sig-boxplots") +
    theme_msig_boxplot()


#### ---- D. PAAD ---- ####
# Associations between KRAS alleles and mutational sigs. in PAAD
# original script: "src/50_35_mutational-signature-allele-associations.R"

panel_D1 <- read_fig_proto("PAAD_median_allele-sig-associate.svg") +
    theme_msig_heatmap() +
    labs(tag = "d")

panel_D2 <- read_fig_proto("PAAD_allele-sig-boxplots") +
    theme_msig_boxplot()


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (
        ((panel_A1 + panel_A2) + plot_layout(widths = c(2, 3))) /
        ((panel_B1 + panel_B2) + plot_layout(widths = c(2, 3))) /
        ((panel_C1 + panel_C2) + plot_layout(widths = c(2, 3))) /
        ((panel_D1 + panel_D2) + plot_layout(widths = c(2, 3)))
    ) +
        plot_layout(heights = c(2, 2, 1, 2))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
