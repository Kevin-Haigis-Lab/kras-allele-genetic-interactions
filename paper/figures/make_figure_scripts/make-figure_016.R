# Figure 016. Supplemental figure to Fig 4 for comutation of LUAD.

FIGNUM <- 16

#> SET THE FIGURE DIMENSIONS
FIG_DIMENSIONS <- get_figure_dimensions(2, "tall")


theme_fig16 <- function(tag_margin = margin(-1, -1, -1, -1, "mm")) {
    theme_comutation() %+replace%
    theme(
        legend.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin)
    )
}

#' Special theme for graphs from 'ggraph'.
theme_graph_fig16 <- function(tag_margin = margin(0, 0, 0, 0, "mm")) {
    theme_graph(base_size = 7, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = tag_margin),
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


#### ---- A. PPI of focal adhesions ---- ####
# The PPI of focal adhesions - an enriched function in the comutation network of
# G12D in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_A <- read_fig_proto("LUAD_G12C_Transcription-Factor-PPIs_MYC") +
    theme_graph_fig16() +
    scale_color_manual(
        drop = FALSE,
        values = c(comut_updown_pal,
                   "none" = "grey70",
                   "in_geneset" = "grey40"),
        label = function(x) { str_replace_all(x, "_", " ") },
        guide = guide_legend(
            title = "comutation interaction",
            nrow = 2,
            title.position = "top",
            title.hjust = 0,
            label.position = "right",
            label.hjust = 0
        )
    ) +
    theme(
        legend.position = c(0.15, 0.95),
        legend.key.height = unit(3, "mm")
    ) +
    labs(tag = "a")



#### ---- B. PPI of Myc ---- ####
# The PPI of Myc (TF) - an enriched function in the comutation network of
# G12C in LUAD.
# original script: "src/20_47_enriched-functions_signaling-pathways.R"

panel_B <- read_fig_proto("LUAD_G12D_KEGG-2019-Human_Focal-adhesion") +
    theme_graph_fig16() +
    theme(
        legend.position = "none"
    ) +
    labs(tag = "b")


#### ---- C. LUAD: KRAS WT vs mut survival curve and table ---- ####
# LUAD: KRAS WT vs mut survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_C_leg_regex <- "kras_allele_grp="
panel_C_pal <- c(WT = "grey50", mut = "black")
names(panel_C_pal) <- paste0(panel_C_leg_regex, names(panel_C_pal))
panel_C <- read_fig_proto("krasmut_survival_LUAD.rds")
panel_C <- fig16_prepare_survplot(panel_C, panel_C_leg_regex, panel_C_pal,
                                   tag = "c",
                                   pval = 0.290,
                                   test_lbl = "log-rank test")


#### ---- D. LUAD: KRAS allele survival curve and table ---- ####
# LUAD: KRAS allele survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_D_leg_regex <- "kras_allele_grp="
panel_D_pal <- short_allele_pal
names(panel_D_pal) <- paste0(panel_D_leg_regex, names(panel_D_pal))
panel_D <- read_fig_proto("krasallele_survival_LUAD.rds")
panel_D <- fig16_prepare_survplot(panel_D, panel_D_leg_regex, panel_D_pal,
                                   tag = "d",
                                   pval = 0.453,
                                   test_lbl = "likelihood ratio test")


#### ---- E. LUAD: KRAS WT vs G12C survival curve and table ---- ####
# LUAD: KRAS WT vs G12C survival curve and table
# original script: "src/70_10_survival-analysis.R"

panel_E_leg_regex <- "kras_allele_grp="
panel_E_pal <- short_allele_pal[c("WT", "G12C")]
names(panel_E_pal) <- paste0(panel_E_leg_regex, names(panel_E_pal))
panel_E <- read_fig_proto("G12C-vs-WT_LUAD.rds")
panel_E <- fig16_prepare_survplot(panel_E, panel_E_leg_regex, panel_E_pal,
                                  tag = "e",
                                  pval = 0.052,
                                  test_lbl = "log-rank test")


#### ---- F-I. Survival curves of comutated genes ---- ####
# Survival curves of genes comutating with G12C.
# original script: "src/70_15_comutation-survival-analysis.R"

proto_paths <- c("survival_alleleorwt_CHRNB4-G12C-LUAD",
                 "survival_alleleorwt_VN1R2-G12C-LUAD",
                 "survival_alleleorwt_ZNF445-G12C-LUAD",
                 "survival_alleleorwt_ZNF804A-G12C-LUAD")


parse_survival_info <- function(cancer, allele, gene) {
    survival_analysis_hits %>%
        filter(cancer == !!cancer &
               interaction_allele == !!allele &
               hugo_symbol == !!gene) %>%
        mutate(
            lbl = paste0("likelihood ratio test: ",
                         str_round(likelihood_ratio_test_pval, 3),
                        "<br>*KRAS* G12C: ",
                        str_round(allele_pval, 3),
                        "<br>*", hugo_symbol, "*: ",
                        str_round(comutation_pval, 3))
        ) %>%
        mutate(x = 240, y = 0.95)
}


prepare_survival_curve <- function(path, title_col = "black") {
    info <- basename(path) %>%
        file_sans_ext() %>%
        str_remove("survival_alleleorwt_") %>%
        str_split_fixed("-", 3)

    info <- parse_survival_info(info[, 3], info[, 2], info[, 1])

    read_fig_proto(path) +
        geom_richtext(
            aes(x = x, y = y, label = lbl),
            data = info, hjust = 1, vjust = 1.0,
            size = 1.8, family = "arial",
            fill = NA, label.color = NA, # remove background and outline
            label.padding = unit(rep(0, 4), "pt") # remove padding
        ) +
        theme_fig16(margin(-1, -1, -1, -2, "mm")) +
        theme(
            plot.margin = margin(0, 2, 0, 0, "mm"),
            plot.title = element_text(size = 7, hjust = 0.5, color = title_col),
            axis.text = element_text(family = "arial", size = 5),
            legend.position = "none",
        ) +
        labs(x = "time (days)", y = "")
}

title_colors <- c(rep(comut_updown_pal[["increased"]], 3),
                  darken(comut_updown_pal[["reduced"]]))
survival_curves <- purrr::map2(proto_paths,
                               title_colors,
                               prepare_survival_curve)

tags <- c("f", "g", "h", "i")
for (i in seq(1, 4)) {
    survival_curves[[i]] <- survival_curves[[i]] + labs(tag = tags[[i]])
}

survival_curves[[1]] <- survival_curves[[1]] + labs(y = "survival probability")

panel_FI <- wrap_plots(survival_curves, nrow = 1)

panel_FI_legend <- read_fig_proto("custom_survival_curve_legend") +
    scale_x_continuous(limits = c(0.7, 4.3)) +
    theme_void() +
    theme(
        legend.position = "none"
    )



#### ---- Figure assembly ---- ####


{
    set.seed(0)

    panel_FI_full <- wrap_elements(
        full = (panel_FI) /
        (
            plot_spacer() + panel_FI_legend + plot_spacer() +
            plot_layout(widths = c(1, 2, 1))
        ) +
        plot_layout(heights = c(10, 1))
    )

    # COMPLETE FIGURE
    full_figure <- ((panel_A | panel_B) + plot_layout(widths = c(2, 1))) /
        (panel_C | panel_D | panel_E) /
        panel_FI_full +
        plot_layout(heights = c(5, 5, 3))

    save_figure(
        full_figure,
        figure_num = FIGNUM,
        dim = FIG_DIMENSIONS
    )
}
