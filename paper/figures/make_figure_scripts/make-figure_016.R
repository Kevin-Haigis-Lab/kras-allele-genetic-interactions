# Figure 016. #> BRIEF DESCRIPTION OF THE FIGURE.

FIGNUM <- 16

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig16 <- function() {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )
}

#' Special theme for graphs from 'ggraph'.
theme_graph_fig16 <- function() {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(0, 0, 0, 0, "mm")),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "bottom",
        legend.title = element_text(size = 5, hjust = 0.5),
        legend.text = element_text(size = 5, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "mm")
    )
}



#### ---- Prepare Survival plots ---- ####


theme_classic_fig16 <- function() {
    theme_classic_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-1, -1, -1, -1, "mm"))
    )

}


# Use 'patchwork' to put the curve and table together.
fig16_patch_ggsurvplot <- function(ggsurv_obj, layout_heights = c(9, 3, 1)) {
    survival_curve <- ggsurv_obj$plot
    survival_tbl <- ggsurv_obj$table
    p <- (survival_curve / survival_tbl / guide_area()) +
        plot_layout(heights = layout_heights, guides = "collect")
    return(p)
}


# Prepare the survplot survival curve and table.
fig16_modify_survplot_components <- function(p, leg_regex, pal,
                                              legend_nrow = 1,
                                              tag = NULL) {
    p$plot <- p$plot +
        scale_color_manual(
            labels = function(x) str_remove(x, leg_regex),
            values = pal,
            guide = guide_legend(title = NULL,
                                 nrow = legend_nrow,
                                 label.position = "top",
                                 label.vjust = -5)
        ) +
        scale_x_continuous(expand = expansion(mult = c(c(0.04, 0.01)))) +
        theme_fig16() +
        theme(
            legend.position = "bottom",
            axis.title.x = element_blank(),
            plot.title = element_blank(),
            legend.box.margin = margin(0, 0, 0, 0, "mm"),
            legend.box.spacing = unit(0, "mm")
        ) +
        labs(
            y = "survival probability",
            tag = tag
        )

    p$table <- p$table +
        scale_color_manual(
            labels = function(x) str_remove(x, leg_regex),
            values = pal
        ) +
        scale_y_discrete(
            labels = function(x) fct_inorder(rev(str_remove(x, leg_regex)))
        ) +
        scale_x_continuous(expand = expansion(mult = c(c(0.04, 0.01)))) +
        theme_classic_fig16() +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            plot.title = element_blank(),
            plot.margin = margin(0, 0, 2, 0, "mm")
        ) +
        labs(x = "time (days)")
    return(p)
}


add_pval_survplot <- function(p, pval, test_lbl) {
    pval <- glue("{test_lbl} p-value: {str_round(pval, 3, 3)}")
    dat <- tibble(label = pval, x = 240, y = 0.95,)
    p$plot <- p$plot +
            geom_richtext(aes(label = label, x = x, y = y),
                          data = dat,
                          hjust = 1, vjust = 1.0,
                          size = 1.8, family = "arial",
                          fill = NA, label.color = NA,
                          label.padding = unit(rep(0, 4), "pt"))
    return(p)
}


fig16_prepare_survplot <- function(p, leg_regex, pal, legend_nrow = 1,
                                   tag = NULL,
                                   pval = NULL, test_lbl = NULL) {
    p <- style_ggsurvminer_plot(p)
    p <- fig16_modify_survplot_components(p, leg_regex, pal, legend_nrow,
                                           tag = tag)
    if (!is.null(pval) & !is.null(test_lbl)) {
        p <- add_pval_survplot(p, pval, test_lbl)
    }

    p <- fig16_patch_ggsurvplot(p)
}


#### ---- A. LUAD: KRAS WT vs mut survival curve and table ---- ####
# LUAD: KRAS WT vs mut survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_A_leg_regex <- "kras_allele_grp="
panel_A_pal <- c(WT = "grey50", mut = "black")
names(panel_A_pal) <- paste0(panel_A_leg_regex, names(panel_A_pal))

panel_A <- read_fig_proto("krasmut_survival_LUAD.rds")
panel_A <- fig16_prepare_survplot(panel_A, panel_A_leg_regex, panel_A_pal,
                                   tag = "a",
                                   pval = 0.290, test_lbl = "log-rank test")


#### ---- B. LUAD: KRAS allele survival curve and table ---- ####
# LUAD: KRAS allele survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_B_leg_regex <- "kras_allele_grp="
panel_B_pal <- short_allele_pal
names(panel_B_pal) <- paste0(panel_B_leg_regex, names(panel_B_pal))

panel_B <- read_fig_proto("krasallele_survival_LUAD.rds")
panel_B <- fig16_prepare_survplot(panel_B, panel_B_leg_regex, panel_B_pal,
                                   tag = "b",
                                   pval = 0.453,
                                   test_lbl = "likelihood ratio test")


#### ---- C. LUAD: KRAS WT vs G12C survival curve and table ---- ####
# LUAD: KRAS WT vs G12C survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_C_leg_regex <- "kras_allele_grp="
panel_C_pal <- short_allele_pal[c("WT", "G12C")]
names(panel_C_pal) <- paste0(panel_C_leg_regex, names(panel_C_pal))
panel_C <- read_fig_proto("G12C-vs-WT_LUAD.rds")
panel_C <- fig16_prepare_survplot(panel_C, panel_C_leg_regex, panel_C_pal,
                                   tag = "c",
                                   pval = 0.052, test_lbl = "log-rank test")


#### ---- D. PPI of Myc ---- ####
# The PPI of Myc (TF) - an enriched function in the comutation network of
# G12C in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_D <- read_fig_proto("LUAD_G12C_Transcription-Factor-PPIs_MYC") +
    theme_graph_fig16() +
    labs(tag = "d")


#### ---- E. PPI of focal adhesions ---- ####
# The PPI of focal adhesions - an enriched function in the comutation network of
# G12D in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_E <- read_fig_proto("LUAD_G12D_KEGG-2019-Human_Focal-adhesion") +
    theme_graph_fig16() +
    labs(tag = "e")


#### ---- Figure assembly ---- ####

{
    set.seed(0)

    # COMPLETE FIGURE
    full_figure <- (panel_A | panel_B  | panel_C) / panel_D / panel_E +
        plot_layout(heights = c(1, 1, 1))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
