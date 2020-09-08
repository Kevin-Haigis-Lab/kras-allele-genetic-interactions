
library(tidyverse)

lib_string <- "stats, glue, conflicted, assertr, testthat,
        glmnet, parallel, caret, ggfortify, tidygraph, jhcutils,
        magrittr, ggpubr, ggraph, ggtext, patchwork, ggplot2, broom,
        tidyverse" %>%
    str_split(",") %>%
    unlist() %>%
    str_squish() %>%
    map(function(a) { library(a, character.only = TRUE) })

for (f in list.files("lib", pattern = "R$", full.names = TRUE)) {
    if (str_detect(f, "enrich|global")) { next }
    source(f)
}

#### ---- Conflicts ---- ####
# Declare which namespaces to use for conflicting functions.

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("cache", "ProjectTemplate")
conflict_prefer("rename", "dplyr")
conflict_prefer("parLapply", "parallel")
conflict_prefer("which", "Matrix")


#### ---- Options ---- ####
options(dplyr.summarise.inform = FALSE)


theme_set(theme_bw(base_size = 11, base_family = "Arial"))

cancer_muts_adjvaf <- readRDS("~/Downloads/cancer_muts_adjvaf.rds")
GRAPHS_DIR <- "90_03_mutation-burden-distribution"


density_margins_plot <- function(df, x, y, x_label, y_label) {

    margin_plot <- function(df, x) {
        df %>%
            ggplot(aes({{ x }})) +
            geom_density(color = "dodgerblue", fill = "dodgerblue",
                         size = 1.2, alpha = 0.4) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
    }

    p1 <- df %>%
        ggplot(aes({{ x }}, {{ y }})) +
        geom_jitter(alpha = 0.3, size = 0.1, width = 0.02, height = 0.02) +
        geom_density2d(color = "tomato", alpha = 0.9, size = 1.2) +
        labs(x = x_label,
             y = y_label)

    p2 <- margin_plot(df, {{ x }})+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
    p3 <- margin_plot(df, {{ y }}) +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank()) +
        coord_flip()


    p_design <- c(
        area(t = 1, l = 1, b = 2, r = 5),
        area(t = 3, l = 1, b = 7, r = 5),
        area(t = 3, l = 6, b = 7, r = 7)
    )
    patch <- p2 + p1 + p3 +
        plot_layout(design = p_design) &
        theme(axis.ticks = element_blank())
    return(patch)
}

density_margins_per_cancer <- function(cancer, data, x, y, x_label, y_label,
                                       plt_name_glue) {
    p <- density_margins_plot(
        data,
        x = {{ x }},
        y = {{ y }},
        x_label = x_label,
        y_label = y_label
    ) +
        plot_annotation(title = cancer,
                        theme = theme(plot.title = element_text(hjust = 0.5)))

    plt_name <- as.character(glue(plt_name_glue))
    ggsave_wrapper(p,
                   plot_path(GRAPHS_DIR, plt_name),
                   "medium")
}

cancer_muts_adjvaf %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(
        density_margins_per_cancer,
        x = purity,
        y = VAF,
        x_label = "tumor sample purity",
        y_label = "VAF (unadjusted)",
        plt_name_glue = "vaf-density2d_{cancer}.svg"
    )

cancer_muts_adjvaf %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(
        density_margins_per_cancer,
        x = purity,
        y = adjusted_vaf,
        x_label = "tumor sample purity",
        y_label = "VAF (adjusted for tumor purity)",
        plt_name_glue = "vaf-adjusted-density2d_{cancer}.svg"
    )




subtitle <- "The dashed line represents where the adjustment for tumor purity had no effect.
Each point represents the mean value for a single tumor sample."

vaf_v_adjvaf_scatter <- cancer_muts_adjvaf %>%
    group_by(cancer, tumor_sample_barcode) %>%
    summarise(avg_VAF = mean(VAF),
              avg_adjusted_vaf = mean(adjusted_vaf),
              purity = mean(purity)) %>%
    ggplot(aes(avg_VAF, avg_adjusted_vaf)) +
    geom_abline(slope = 1, intercept = 0, lty = 2, color = "grey50") +
    geom_jitter(aes(size = purity, color = cancer),
                alpha = 0.5, height = 0.02, width = 0.02) +
    scale_size_continuous(range = c(1, 5)) +
    scale_color_manual(values = cancer_palette, drop = TRUE) +
    labs(x = "VAF",
         y = "VAF adjusted for tumor sample purity",
         size = "tumor sample purity",
         color = "cancer type",
         title = "Effects of adjusting VAF for tumor sample purity",
         subtitle = subtitle)
ggsave_wrapper(vaf_v_adjvaf_scatter,
               plot_path(GRAPHS_DIR, "vaf-v-adjvaf-scatter.svg"),
               "medium")


jitter_pos <- position_jitter(width = 0.02, height = 0, seed = 0)
subtitle <- "Each point represents the median (± the std. dev.) VAF of an individual tumor sample."
purity_vs_vaf_scatter <- cancer_muts_adjvaf %>%
    group_by(cancer, tumor_sample_barcode) %>%
    summarise(avg_VAF = median(VAF),
              sd_VAF = sd(VAF),
              purity = mean(purity)) %>%
    ungroup() %>%
    mutate(upper_VAF = minmax(avg_VAF + sd_VAF, 0, 1),
           lower_VAF = minmax(avg_VAF - sd_VAF, 0, 1)) %>%
    ggplot(aes(x = purity, y = avg_VAF)) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_linerange(aes(ymin = lower_VAF, ymax = upper_VAF),
                   alpha = 0.3, position = jitter_pos) +
    geom_point(position = jitter_pos) +
    geom_smooth(method = "lm", formula = "y ~ x") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    labs(x = "tumor sample purity",
         y = "mutation VAF (unadjusted for tumor sample purity)",
         title = "Correlation between tumor sample purity and VAF values",
         subtitle = subtitle)
ggsave_wrapper(purity_vs_vaf_scatter,
               plot_path(GRAPHS_DIR, "purity-vs-vaf-scatter.svg"),
               "large")





# TODO: compare VAF values with synonymous mutations
# Particular focus on KRAS and its comutation partners


cancer_silent_muts_df <- readRDS("~/Downloads/cancer_silent_muts_df.rds")
