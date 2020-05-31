# Figure 045. #> The first group of results from the analysis for explaining
# allele-specific genetic dependency by comutation.

FIGNUM <- 45

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(1, "medium")


theme_fig45 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


theme_coefplot <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
    theme_fig45(tag_margin = tag_margin) %+replace%
    theme(
        plot.title = element_markdown(size = 7, face = "bold", family = "Arial"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown()
    )
}


theme_lineplot <- function(tag_margin = margin(1, 1, 1, 1, "mm")) {
    theme_minimal_comutation() %+replace%
    theme(
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        panel.border = element_blank()
    )
}


theme_sideline <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_void(base_size = 7, base_family = "Arial") %+replace%
    theme(
        axis.title.y = element_markdown(angle = 90),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}


#### ---- A. TP53 explaining STARD9 ---- ####
# TP53 explaining dep. of G12D on STARD9 in COAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_A_1 <- read_fig_proto("COAD_G12D_STARD9_coef-plot") +
    theme_coefplot()
panel_A_2 <- read_fig_proto("COAD_G12D_STARD9_line-plot") +
    theme_lineplot()

panel_A_sideline <- tibble(x = 1, y = c(1, 2)) %>%
    ggplot(aes(x, y)) +
    geom_line(aes(group = x), size = 0.7, color = "grey50",
              lineend = "round") +
    scale_y_continuous(limits = c(1, 2), expand = c(0.01, 0.01)) +
    theme_sideline() +
    labs(y = "Effect of *TP53* mutation on dependency on *STARD9* in COAD",
         tag = "a")



#### ---- B. SMAD4 explaining EEF1E1 ---- ####
# SMAD4 explaining dep. of G12D on EEF1E1 in PAAD.
# original script: "src/40_63_explore-synlet-comut.R"

panel_B_1 <- read_fig_proto("PAAD_G12D_EEF1E1_coef-plot") +
    theme_coefplot()
panel_B_2 <- read_fig_proto("PAAD_G12D_EEF1E1_line-plot") +
    theme_lineplot()

panel_B_sideline <- tibble(x = 1, y = c(1, 2)) %>%
    ggplot(aes(x, y)) +
    geom_line(aes(group = x), size = 0.7, color = "grey50",
              lineend = "round") +
    scale_y_continuous(limits = c(1, 2), expand = c(0.01, 0.01)) +
    theme_sideline() +
    labs(y = "Effect of *SMAD4* mutation on dependency on *EEF1E1* in PAAD",
         tag = "b")


#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (
        ((panel_A_sideline | (panel_A_1 / panel_A_2)) +
            plot_layout(widths = c(1, 20))) /
        plot_spacer() /
        ((panel_B_sideline | (panel_B_1 / panel_B_2)) +
            plot_layout(widths = c(1, 20)))
    ) +
        plot_layout(heights = c(20, 1, 20))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
