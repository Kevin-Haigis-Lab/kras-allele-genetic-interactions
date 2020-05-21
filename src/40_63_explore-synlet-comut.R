# Explore some of the results from the modeling of synthetic lethality
# using comutation interactions.

GRAPHS_DIR <- "40_63_explore-synlet-comut"
# reset_graph_directory(GRAPHS_DIR)

################################################################################
## REMOVE FOR O2 ##

library(stats)
library(glue)
library(conflicted)
library(ggfortify)
library(tidygraph)
library(jhcutils)
library(magrittr)
library(ggpubr)
library(ggraph)
library(ggtext)
library(patchwork)
library(ggplot2)
library(broom)
library(tidyverse)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("cache", "ProjectTemplate")
conflict_prefer("rename", "dplyr")
conflict_prefer("parLapply", "parallel")
conflict_prefer("which", "Matrix")


load("cache/synlet_comut_model_res.RData")

synlet_comut_model_res %<>%
    filter(
        (cancer == "COAD" & allele == "G12D" & hugo_symbol == "STARD9") |
            (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "EEF1E1") |
            (cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5") |
            (cancer == "PAAD" & allele == "G12R" & hugo_symbol == "KIAA1257") |
            (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "DMBX1")
    )


source("lib/ggplot2_helpers.R")
source("lib/helpers.R")

################################################################################



pastel_red <- "#FF8FAC"
pastel_blue <- "#4096B3"



#### ---- MASKING ---- ####

# Make the coefficient plot for the masking group.
make_mask_coef_plot <- function(cancer, allele, hugo_symbol, fit, ...) {
    p <- fit$elastic_model %>%
        broom::tidy() %>%
        filter(term != "(Intercept)") %>%
        mutate(term = ifelse(term == "kras_allele", !!allele, term),
               term = fct_reorder(term, estimate)) %>%
        ggplot(aes(x = estimate, y = term)) +
        geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
        geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
        geom_label(aes(label = term, fill = estimate),
                   family = "Arial", size = 2.2,
                   label.padding = unit(0.8, "mm"),
                   label.r = unit(0.4, "mm"),
                   color = "black", label.size = 0) +
        scale_color_gradient2(high = pastel_red, mid = "grey90", low = pastel_blue) +
        scale_fill_gradient2(high = pastel_red, mid = "grey90", low = pastel_blue) +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            plot.title = element_blank(),
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_markdown()
        ) +
        labs(x = glue(
            "estimated effect of mutation on dependency on *{hugo_symbol}*"
        ))
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR,
                  glue("{cancer}_{allele}_{hugo_symbol}_coef-plot.svg")),
        "small"
    )
    return(p)
}


mask_pal <- c(short_allele_pal["G12D"],
              TP53 = "#FF9945",
              SMAD4 = "#AACC7C")

# Make the line plot of gene effect for the masking group.
make_mask_line_plot <- function(cancer, allele, hugo_symbol, fit, other_gene, ...) {
    p_data <- fit$data %>%
        rename(other_gene_col = !!other_gene) %>%
        filter(xor(kras_allele == 1, other_gene_col == 1)) %>%
        mutate(mutation = case_when(
            kras_allele == 1 ~ "G12D",
            other_gene_col == 1 ~ !!other_gene
        )) %>%
        mutate(mutation = fct_rev(factor(mutation)))

    label_data <- p_data %>%
        group_by(mutation) %>%
        summarise(x = mean(gene_effect)) %>%
        ungroup() %>%
        mutate(y = 1.7)

    p <- p_data %>%
        mutate(y = ifelse(mutation == !!allele, 0.8, 1.2)) %>%
        ggplot(aes(x = gene_effect, y = y, color = mutation)) +
        geom_hline(yintercept = c(0.8, 1.2), size = 0.1, color = "grey70") +
        geom_point(size = 5, alpha = 0.7) +
        geom_richtext(data = label_data,
                      aes(x = x, y = y, label = mutation),
                      size = 3, family = "Arial", fontface = "bold",
                      fill = NA, label.color = NA,
                      label.padding = unit(rep(0, 4), "mm")) +
        scale_color_manual(values = mask_pal, drop = TRUE) +
        scale_y_continuous(limits = c(0.5, 2.1),
                           expand = expansion(mult = c(0, 0))) +
        theme_minimal(base_size = 7, base_family = "Arial") +
        theme(
            legend.position = "none",
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_markdown()
        ) +
        labs(x = glue("*{hugo_symbol}* dependency score"),
             y = NULL)

    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR,
                  glue("{cancer}_{allele}_{hugo_symbol}_line-plot.svg")),
        "small"
    )
    return(p)
}


# Make the patchwork for the masking group.
make_mask_patch_plot <- function(cancer, allele, hugo_symbol, fit, other_gene, ...) {
    coef_p <- make_mask_coef_plot(cancer, allele, hugo_symbol, fit)
    line_p <- make_mask_line_plot(cancer, allele, hugo_symbol, fit, other_gene)
    patch <- (coef_p / line_p) +
        plot_layout(heights = c(2, 1))
    ggsave_wrapper(
        patch,
        plot_path(GRAPHS_DIR, glue("{cancer}_{allele}_{hugo_symbol}_patch.svg")),
        "small"
    )
}


# All of the groups to make masking group plots for.
masking_hits <- tribble(
    ~ cancer, ~ allele, ~ hugo_symbol, ~ other_gene,
    "COAD", "G12D", "STARD9", "TP53",
    "PAAD", "G12D", "EEF1E1", "SMAD4",
    "PAAD", "G12D", "MYBL2", "SMAD4",
    "PAAD", "G12D", "ABI1", "SMAD4",
)

synlet_comut_model_res %>%
    right_join(masking_hits, by = c("cancer", "allele", "hugo_symbol")) %>%
    pwalk(make_mask_patch_plot)









################################################################################





#### ---- COAD G12D: SRSF5 ---- ####

# G12D has reduced comut with HECW1 and increased comut with APC.
# When all three are mutated, there is less dep. on SRSF5.

# group: linking


srsf5_res <- synlet_comut_model_res %>%
    filter(cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5")

srsf5_comuts <- c("SORL1", "APC", "HECW1")
srsf5_terms <- c("kras_allele", srsf5_comuts)
srsf5_coef_plot <- srsf5_res$fit[[1]]$elastic_model %>%
    broom::tidy(return_zeros = TRUE) %>%
    filter(term %in% srsf5_terms) %>%
    mutate(term = ifelse(term == "kras_allele", !!srsf5_res$allele, term),
           term = fct_reorder(term, estimate)) %>%
    ggplot(aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
    geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
    geom_label(aes(label = term, fill = estimate),
               family = "Arial", size = 2.2,
               label.padding = unit(0.8, "mm"),
               label.r = unit(0.4, "mm"),
               color = "black", label.size = 0) +
    scale_color_gradient2(high = pastel_red, mid = "grey90", low = pastel_blue) +
    scale_fill_gradient2(high = pastel_red, mid = "grey90", low = pastel_blue) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        panel.border = element_blank()
    ) +
    labs(x = glue(
        "estimated effect of mutation on dependency on *{srsf5_res$hugo_symbol[[1]]}*"
    ))


srsf5_box_data <- tibble()
for (gene in srsf5_comuts) {
    srsf5_box_data <- bind_rows(
        srsf5_res$fit[[1]]$data %>%
            add_column(comut_gene = tidyselect::all_of(!!gene)) %>%
            rename(highlight = gene) %>%
            select(gene_effect, kras_allele, comut_gene, highlight),
        srsf5_box_data
    )
}

srsf5_rects <- tibble(
    xmin = -Inf,
    xmax = Inf,
    ymin = c(0.2, 0),
    ymax = c(0, -0.61),
    fill = c("skyblue", "grey90")
)

srsf5_box_plot <- srsf5_box_data %>%
    mutate(kras_allele = ifelse(kras_allele == 1, !!srsf5_res$allele, "other"),
           highlight = ifelse(highlight == 1, "mutant", "WT")) %>%
    ggplot() +
    geom_rect(data = srsf5_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = fill),
              alpha = 0.2) +
    geom_hline(yintercept = 0, size = 0.2, color = "grey30") +
    geom_jitter(aes(x = comut_gene, y = gene_effect,
                        color = highlight, shape = kras_allele), width = 0.2) +
    scale_color_manual(values = c("grey20", "grey75")) +
    scale_shape_manual(values = c(G12D = 17, other = 16)) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(-0.61, 0.2), expand = c(0, 0)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.ticks = element_blank(),
        legend.title = element_markdown()
    ) +
    labs(y = glue("*{srsf5_res$hugo_symbol}* dependency score"),
         color = "comutation",
         shape = "*KRAS* allele")


srsf5_comut_a <- tibble(
    y = 1,
    x = 2,
    label = c("G12D")
)
srsf5_comut_b <- tibble(
    x = c(2, 1, 1, 2, 2, 2, 2, 3, 3),
    y = rep(c(1, 1.5, 3), 3),
    group = rep(srsf5_comuts, each = 3),
    interaction = rep(c("increased", "reduced", "reduced"), each = 3)
)
srsf5_comut_c <- tibble(
    x = c(1, 3),
    y = c(1.5, 1.5),
    label = c("increased", "reduced")
)

srsf5_comut_plot <- ggplot() +
    ggforce::geom_bezier(
        aes(x = x, y = y, group = group, color = interaction),
        data = srsf5_comut_b
    ) +
    geom_label(
        aes(x = x, y = y, label = label),
        data = srsf5_comut_a,
        family = "Arial", size = 2.2,
        label.padding = unit(0.8, "mm"),
        label.r = unit(0.4, "mm"),
        color = "black", label.size = 0
    ) +
    geom_text(
        aes(x = x, y = y, label = label, color = label),
        data = srsf5_comut_c,
        family = "Arial", size = 2.2,
    ) +
    scale_x_continuous(limits = c(0.5, 3.5)) +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0))) +
    scale_color_manual(values = comut_updown_pal,
                       guide = FALSE) +
    theme_void(base_size = 7, base_family = "Arial")


srsf5_patch <- (srsf5_coef_plot / srsf5_box_plot / srsf5_comut_plot) +
    plot_layout(heights = c(3, 9, 1))
ggsave_wrapper(
    srsf5_patch,
    plot_path(GRAPHS_DIR, "COAD_G12D_SRSF5-patch.svg"),
    "small"
)









#### ---- PAAD G12R: KIAA1257 ---- ####

# G12D has reduced comut with DNAH5.
# When one is mutated, dep. on KIAA1257 is increased.
# When both are mutated, dep. on KIAA1257 is even stronger.

# group: linking

kiaa1257_res <- synlet_comut_model_res %>%
    filter(cancer == "PAAD" & allele == "G12R" & hugo_symbol == "KIAA1257")

kiaa1257_pal <- c("neither" = "grey50",
                  "DNAH5" = "#F59237",
                  "G12R" = short_allele_pal[["G12R"]],
                  "G12R & DNAH5" = "#FF6B72")

p_data <- kiaa1257_res$fit[[1]]$data %>%
    mutate(mutation = case_when(
        kras_allele + DNAH5 == 2 ~ "G12R & DNAH5",
        kras_allele == 1 ~ "G12R",
        DNAH5 == 1 ~ "DNAH5",
        TRUE ~ "neither"
    ),
    mutation = factor(mutation, names(kiaa1257_pal)))

p_data_summary <- p_data %>%
    group_by(mutation) %>%
    summarise(gene_effect = mean(gene_effect)) %>%
    ungroup()

set.seed(123)
kiaa1257_box_plot <- p_data %>%
    ggplot(aes(y = gene_effect, x = mutation)) +
    geom_line(data = p_data_summary,
              group = 1, color = "grey50", alpha = 0.5) +
    geom_point(data = p_data_summary,
               aes(color = mutation),
               size = 3, shape = 18) +
    geom_jitter(aes(color = mutation), height = 0, width = 0.2, alpha = 0.8) +
    scale_color_manual(values = kiaa1257_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown()
    ) +
    labs(y = glue(
        "dependency score of *{kiaa1257_res$hugo_symbol[[1]]}*"
    ))
ggsave_wrapper(
    kiaa1257_box_plot,
    plot_path(GRAPHS_DIR, "PAAD_G12R_KIAA1257_box-plot.svg"),
    "small"
)


kiaa1257_coef_plot <- kiaa1257_res$fit[[1]]$elastic_model %>%
    broom::tidy() %>%
    filter(term != "(Intercept)") %>%
    mutate(term = ifelse(term == "kras_allele", kiaa1257_res$fit[[1]]$allele, term)) %>%
    ggplot(aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
    geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
    geom_label(aes(label = term, fill = estimate),
               family = "Arial", size = 2.2,
               label.padding = unit(0.8, "mm"),
               label.r = unit(0.4, "mm"),
               color = "black", label.size = 0) +
    scale_color_gradient2(high = pastel_red, mid = "grey90", low = pastel_blue) +
    scale_fill_gradient2(high = pastel_red, mid = "grey90", low = pastel_blue) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        panel.border = element_blank()
    ) +
    labs(x = glue(
        "estimated effect of mutation on dependency on *{kiaa1257_res$hugo_symbol[[1]]}*"
    ))
ggsave_wrapper(
    kiaa1257_coef_plot,
    plot_path(GRAPHS_DIR, "PAAD_G12R_KIAA1257_coef-plot.svg"),
    "small"
)


g12r_dnah5_line <- tibble(name = c("DNAH5", "reduced comutation", "G12R"),
                          color = c("white", comut_updown_pal[["reduced"]], "white"),
                          y = c(1, 1.3, 1)) %>%
    mutate(name = fct_inorder(name)) %>%
    ggplot(aes(x = name, y = y)) +
    geom_segment(x = "DNAH5", xend = "G12R", y = 1, yend = 1,
                 color = comut_updown_pal[["reduced"]],
                 size = 1.4) +
    geom_label(aes(label = name, fill = name, color = color),
               family = "Arial", size = 2.2,
               label.size = 0) +
    scale_color_identity() +
    scale_fill_manual(values = kiaa1257_pal, guide = FALSE) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_void(base_size = 7, base_family = "Arial")


kiaa1257_layout <- c(
    area(t = 1, l = 1, b = 4, r = 10),
    area(t = 5, l = 1, b = 14, r = 10),
    area(t = 7, l = 5.5, b = 8, r = 10)
)

kiaa1257_patch <- kiaa1257_coef_plot + kiaa1257_box_plot + g12r_dnah5_line +
    plot_layout(design = kiaa1257_layout)
ggsave_wrapper(
    kiaa1257_patch,
    plot_path(GRAPHS_DIR, "PAAD_G12R_KIAA1257_patch.svg"),
    "small"
)





#### ---- PAAD G12D: DMBX1 ---- ####

# G12D has increased comut with GPR98 and RNF43.
# When one is mutated, dep. on DMBX1 is increased.
# When KRAS and one other are mutated, dep. on DMBX1 is even stronger.

# group: linking

dmbx1_res <- synlet_comut_model_res %>%
    filter(cancer == "PAAD" & allele == "G12D" & hugo_symbol == "DMBX1")

dmbx1_comuts <- c("GPR98", "RNF43")
