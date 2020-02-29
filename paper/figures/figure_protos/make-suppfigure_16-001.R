# Supplemental Figure 16.  Survival curves of LUAD patients
# to accompany those in Fig 3.

FIGNUM <- 16
SUPPLEMENTAL <- TRUE
VERSION <- 1

FIG_DIMENSIONS <- get_figure_dimensions(2, "short")
FIG_DIMENSIONS$height <- 0.5 * FIG_DIMENSIONS$height


theme_figS16 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}


theme_classic_figS16 <- function() {
    theme_classic_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )

}


# Use 'patchwork' to put the curve and table together.
figS16_patch_ggsurvplot <- function(ggsurv_obj, layout_heights = c(9, 3, 1)) {
    survival_curve <- ggsurv_obj$plot
    survival_tbl <- ggsurv_obj$table
    p <- (survival_curve / survival_tbl / guide_area()) +
        plot_layout(heights = layout_heights, guides = "collect")
    return(p)
}


# Prepare the survplot survival curve and table.
figS16_modify_survplot_components <- function(p, leg_regex, pal,
                                              legend_nrow = 1) {
    p$plot <- p$plot +
        scale_color_manual(
            labels = function(x) str_remove(x, leg_regex),
            values = pal,
            guide = guide_legend(title = NULL,
                                 nrow = legend_nrow,
                                 label.position = "top",
                                 label.vjust = -5)
        ) +
        scale_x_continuous(expand = expand_scale(mult = c(c(0.04, 0.01)))) +
        theme_figS16() +
        theme(
            legend.position = "bottom",
            axis.title.x = element_blank(),
            plot.title = element_blank()
        ) +
        labs(y = "survival probability")

    p$table <- p$table +
        scale_color_manual(
            labels = function(x) str_remove(x, leg_regex),
            values = pal
        ) +
        scale_y_discrete(
            labels = function(x) str_remove(x, leg_regex)
        ) +
        scale_x_continuous(expand = expand_scale(mult = c(c(0.04, 0.01)))) +
        theme_classic_figS16() +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            plot.title = element_blank()
        ) +
        labs(x = "time (days)")
    return(p)
}


figS16_prepare_survplot <- function(p, leg_regex, pal, legend_nrow = 1) {
    p <- style_ggsurvminer_plot(p)
    p <- figS16_modify_survplot_components(p, leg_regex, pal, legend_nrow)
    p <- figS16_patch_ggsurvplot(p)
}




#### ---- A. LUAD: KRAS WT vs mut survival curve and table ---- ####
# LUAD: KRAS WT vs mut survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_A_leg_regex <- "kras_allele_grp="
panel_A_pal <- c(WT = "grey50", mut = "black")
names(panel_A_pal) <- paste0(panel_A_leg_regex, names(panel_A_pal))

panel_A <- read_fig_proto("krasmut_survival_LUAD.rds",
                          FIGNUM, supp = SUPPLEMENTAL)
panel_A <- figS16_prepare_survplot(panel_A, panel_A_leg_regex, panel_A_pal)


#### ---- B. LUAD: KRAS allele survival curve and table ---- ####
# LUAD: KRAS allele survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_B_leg_regex <- "kras_allele_grp="
panel_B_pal <- short_allele_pal
names(panel_B_pal) <- paste0(panel_B_leg_regex, names(panel_B_pal))

panel_B <- read_fig_proto("krasallele_survival_LUAD.rds",
                          FIGNUM, supp = SUPPLEMENTAL)
panel_B <- figS16_prepare_survplot(panel_B, panel_B_leg_regex, panel_B_pal)


#### ---- A. LUAD: KRAS WT vs G12C survival curve and table ---- ####
# LUAD: KRAS WT vs G12C survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_C_leg_regex <- "kras_allele_grp="
panel_C_pal <- short_allele_pal[c("WT", "G12C")]
names(panel_C_pal) <- paste0(panel_C_leg_regex, names(panel_C_pal))
panel_C <- read_fig_proto("G12C-vs-WT_LUAD.rds",
                          FIGNUM, supp = SUPPLEMENTAL)
panel_C <- figS16_prepare_survplot(panel_C, panel_C_leg_regex, panel_C_pal)



#### ---- Figure assembly ---- ####

{
    # COMPLETE FIGURE
    full_figure <- (panel_A | panel_B  | panel_C) +
        plot_layout(widths = c(1, 1, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        supp = SUPPLEMENTAL,
        dim = FIG_DIMENSIONS
    )
}
