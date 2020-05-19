# Explore some of the results from the modeling of synthetic lethality
# using comutation interactions.

GRAPHS_DIR <- "40_63_explore-synlet-comut"
reset_graph_directory(GRAPHS_DIR)

################################################################################
## REMOVE FOR O2 ##

# library(stats)
# library(glue)
# library(conflicted)
# library(ggfortify)
# library(tidygraph)
# library(jhcutils)
# library(magrittr)
# library(ggpubr)
# library(ggraph)
# library(ggtext)
# library(patchwork)
# library(ggplot2)
# library(broom)
# library(tidyverse)
#
# conflict_prefer("select", "dplyr")
# conflict_prefer("filter", "dplyr")
# conflict_prefer("slice", "dplyr")
# conflict_prefer("setdiff", "dplyr")
# conflict_prefer("intersect", "dplyr")
# conflict_prefer("cache", "ProjectTemplate")
# conflict_prefer("rename", "dplyr")
# conflict_prefer("parLapply", "parallel")
# conflict_prefer("which", "Matrix")
#
#
# load("cache/synlet_comut_model_res.RData")
#
# synlet_comut_model_res %<>%
#     filter(
#         (cancer == "COAD" & allele == "G12D" & hugo_symbol == "STARD9") |
#             (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "EEF1E1") |
#             (cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5") |
#             (cancer == "PAAD" & allele == "G12R" & hugo_symbol == "KIAA1257")
#     )
#
#
# source("lib/ggplot2_helpers.R")
# source("lib/helpers.R")

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

# group: interaction


srsf5_res <- synlet_comut_model_res %>%
    filter(cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5")


srsf5_res$fit[[1]]$data %>%
    mutate(mutation = case_when(
        kras_allele + HECW1 + APC == 3 ~ "G12D & APC & HECW1",
        kras_allele + HECW1  == 2 ~ "G12D & HECW1",
        kras_allele + APC  == 2 ~ "G12D & APC",
        kras_allele == 1 ~ "G12D",
        HECW1 == 1 ~ "HECW1",
        APC == 1 ~ "APC",
        TRUE ~ "none"
    )) %>%
    ggplot(aes(x = gene_effect, y = mutation)) +
    geom_jitter(height = 0.1, width = 0)





#### ---- PAAD G12R: KIAA1257 ---- ####

# G12D has reduced comut with DNAH5.
# When one is mutated, dep. on KIAA1257 is increased.
# When both are mutated, dep. on KIAA1257 is even stronger.

# group: interaction


kiaa1257_res <- synlet_comut_model_res %>%
    filter(cancer == "PAAD" & allele == "G12R" & hugo_symbol == "KIAA1257")

kiaa1257_res$fit[[1]]$data %>%
    mutate(mutation = case_when(
        kras_allele + DNAH5 == 2 ~ "G12R & DNAH5",
        kras_allele == 1 ~ "G12R",
        DNAH5 == 1 ~ "DNAH5",
        TRUE ~ "none"
    )) %>%
    ggplot(aes(x = gene_effect, y = mutation)) +
    geom_jitter(height = 0.1, width = 0)
