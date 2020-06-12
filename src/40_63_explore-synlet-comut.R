# Explore some of the results from the modeling of synthetic lethality
# using comutation interactions.

GRAPHS_DIR <- "40_63_explore-synlet-comut"
reset_graph_directory(GRAPHS_DIR)

pastel_red <- "#FF8FAC"
pastel_blue <- "#4096B3"

set.seed(0)



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

options(dplyr.summarise.inform = FALSE)

load("cache/synlet_comut_model_res.RData")

synlet_comut_model_res %<>%
    filter(
        (cancer == "COAD" & allele == "G12D" & hugo_symbol == "STARD9") |
            (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "EEF1E1") |
            (cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5") |
            (cancer == "PAAD" & allele == "G12R" & hugo_symbol == "KIAA1257") |
            (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "FKBP1A")
    )


saveFigRds <- function(...) { invisible(NULL) }

source("lib/ggplot2_helpers.R")
source("lib/helpers.R")

################################################################################



#### ---- MASKING ---- ####

# Make the coefficient plot for the masking group.
make_mask_coef_plot <- function(cancer, allele, hugo_symbol, fit, ...) {
    p <- broom::tidy(fit$elastic_model) %>%
        filter(term != "(Intercept)") %>%
        mutate(
            term = str_replace(term, "kras_allele", !!allele),
            term = str_replace(term, ":", " & "),
            term = fct_reorder(term, estimate)
        ) %>%
        ggplot(aes(x = estimate, y = term)) +
        geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
        geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
        geom_label(aes(label = term, fill = estimate),
                   family = "Arial", size = 2.2,
                   label.padding = unit(0.8, "mm"),
                   label.r = unit(0.4, "mm"),
                   color = "black", label.size = 0) +
        scale_color_gradient2(high = pastel_red,
                              mid = "grey90",
                              low = pastel_blue) +
        scale_fill_gradient2(high = pastel_red,
                             mid = "grey90",
                             low = pastel_blue) +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        theme_minimal(base_size = 7, base_family = "Arial") +
        theme(
            plot.title = element_markdown(hjust = 0.5),
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(title = glue(
            "estimated effect of mutation on dependency on *{hugo_symbol}*"
        ))
    plt_name <- glue("{cancer}_{allele}_{hugo_symbol}_coef-plot.svg") %>%
        as.character()
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, plt_name),
        "small"
    )
    saveFigRds(p, plt_name)
    return(p)
}


mask_pal <- c(short_allele_pal["G12D"],
              TP53 = "#FF9945",
              SMAD4 = "#AACC7C")

# Make the line plot of gene effect for the masking group.
make_mask_line_plot <- function(cancer, allele, hugo_symbol, fit, other_gene,
                                ...) {
    p_data <- fit$data %>%
        rename(other_gene_col = !!other_gene) %>%
        filter(xor(kras_allele == 1, other_gene_col == 1)) %>%
        mutate(mutation = case_when(
            kras_allele == 1 ~ "G12D",
            other_gene_col == 1 ~ !!other_gene
        )) %>%
        mutate(mutation = fct_rev(factor(mutation)))

    p <- p_data %>%
        ggplot(aes(x = mutation, y = gene_effect, color = mutation))

    if (0 < max(p_data$gene_effect) & 0 > min(p_data$gene_effect)) {
        p <- p +
            geom_hline(yintercept = 0, size = 0.2, color = "grey30")
    }

    p <- p +
        ggbeeswarm::geom_quasirandom(size = 3, alpha = 0.7, groupOnX = TRUE) +
        scale_color_manual(values = mask_pal, drop = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            legend.position = "none",
            axis.title.y = element_markdown(),
            axis.title.x = element_blank()
        ) +
        labs(y = glue("*{hugo_symbol}* dependency score"))

    plt_name <- glue("{cancer}_{allele}_{hugo_symbol}_line-plot.svg") %>%
        as.character()
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, plt_name),
        "small"
    )
    saveFigRds(p, plt_name)
    return(p)
}


# Make the patchwork for the masking group.
make_mask_patch_plot <- function(cancer, allele, hugo_symbol, fit, other_gene,
                                 ...) {
    coef_p <- make_mask_coef_plot(cancer, allele, hugo_symbol, fit)
    line_p <- make_mask_line_plot(cancer, allele, hugo_symbol, fit, other_gene)
    patch <- (coef_p / line_p) +
        plot_layout(heights = c(1, 3))
    plt_name <- as.character(glue("{cancer}_{allele}_{hugo_symbol}_patch.svg"))
    ggsave_wrapper(
        patch,
        plot_path(GRAPHS_DIR, plt_name),
        "small"
    )
    saveFigRds(patch, plt_name)
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
################################################################################
################################################################################


#### ---- COAD G12D: SRSF5 ---- ####

# G12D has reduced comut with HECW1 and increased comut with APC.
# When all three are mutated, there is less dep. on SRSF5.

# group: linking

srsf5_res <- synlet_comut_model_res %>%
    filter(cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5")

srsf5_comuts <- c("APC", "HECW1")
srsf5_terms <- c("kras_allele", srsf5_comuts, "kras_allele:APC")
srsf5_coef_plot <- srsf5_res$fit[[1]]$elastic_model %>%
    broom::tidy(return_zeros = TRUE) %>%
    filter(term %in% srsf5_terms) %>%
    mutate(term = str_replace(term, "kras_allele", srsf5_res$fit[[1]]$allele),
           term = str_replace(term, ":", " & "),
           term = fct_reorder(term, estimate)) %>%
    ggplot(aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
    geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
    geom_label(aes(label = term, fill = estimate),
               family = "Arial", size = 2.2,
               label.padding = unit(0.8, "mm"),
               label.r = unit(0.4, "mm"),
               color = "black", label.size = 0) +
    scale_color_gradient2(high = pastel_red,
                          mid = "grey90",
                          low = pastel_blue) +
    scale_fill_gradient2(high = pastel_red,
                         mid = "grey90",
                         low = pastel_blue) +
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
        "estimated effect of mutation on dependency on *SRSF5*"
    ))
saveFigRds(srsf5_coef_plot, "COAD_G12D_SRSF5_coef-plot.rds")


srsf_box_pal <- c(short_allele_pal, other = "grey50")

srsf5_box_plot <- srsf5_res$fit[[1]]$data %>%
    mutate(
        kras_allele = ifelse(kras_allele == 1, !!srsf5_res$allele, "other"),
        mutation = case_when(
            APC + HECW1 == 2 ~ "APC & HECW1",
            APC == 1 ~ "APC",
            HECW1 == 1 ~ "HECW1",
            TRUE ~ "WT"),
        mutation = factor(mutation,
                          levels = c("WT", "APC", "HECW1", "APC & HECW1"))
    ) %>%
    ggplot(aes(mutation, gene_effect)) +
    geom_hline(yintercept = 0, size = 0.2, color = "grey30") +
    geom_jitter(aes(color = kras_allele), width = 0.2, alpha = 0.7) +
    scale_color_manual(values = srsf_box_pal, drop = TRUE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.ticks = element_blank(),
        legend.title = element_markdown()
    ) +
    labs(y = glue("*{srsf5_res$hugo_symbol}* dependency score"),
         color = "*KRAS* allele")
saveFigRds(srsf5_box_plot, "COAD_G12D_SRSF5_box-plot2.rds")


load("cache/genetic_interaction_gr.RData")
srsf5_comut_gr <- genetic_interaction_gr %E>%
    filter(cancer == "COAD") %>%
    mutate(genetic_interaction = switch_comut_terms(genetic_interaction),
           edge_label = paste(genetic_interaction, " comutation")) %N>%
    mutate(name = str_remove(name, "KRAS_")) %>%
    filter(name %in% c("G12D", srsf5_comuts)) %>%
    mutate(node_text_color = ifelse(name %in% kras_dark_lbls, "white", "black"))

srsf5_comut_gr_layout <- create_layout(srsf5_comut_gr, layout = "nicely")
srsf5_comut_gr_layout$x <- c(2, 3, 1)
srsf5_comut_gr_layout$y <- rep(0, 3)
srsf5_comut_plot <- ggraph(srsf5_comut_gr_layout) +
    geom_edge_link(aes(color = genetic_interaction,
                       label = edge_label),
                   vjust = -1, label_size = 2.2, width = 1) +
    geom_node_label(aes(label = name, fill = name, color = node_text_color),
                    family = "Arial", size = 2.2,
                    label.padding = unit(0.8, "mm"),
                    label.r = unit(0.4, "mm"),
                    label.size = 0) +
    scale_edge_color_manual(values = comut_updown_pal, guide = FALSE) +
    scale_color_identity(guide = FALSE) +
    scale_fill_manual(values = short_allele_pal, na.value = "grey90", guide = FALSE) +
    theme_graph() +
    theme(
        plot.margin = margin(3, 0, 0, 0, unit = "mm")
    )
saveFigRds(srsf5_comut_plot, "COAD_G12D_SRSF5_comut-graph-plot.rds")

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
    geom_hline(yintercept = 0, size = 0.2, color = "grey30") +
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
saveFigRds(kiaa1257_box_plot, "PAAD_G12R_KIAA1257_box-plot.rds")
ggsave_wrapper(
    kiaa1257_box_plot,
    plot_path(GRAPHS_DIR, "PAAD_G12R_KIAA1257_box-plot.svg"),
    "small"
)


kiaa1257_coef_plot <- kiaa1257_res$fit[[1]]$elastic_model %>%
    broom::tidy() %>%
    filter(term != "(Intercept)") %>%
    mutate(term = ifelse(term == "kras_allele",
                         kiaa1257_res$fit[[1]]$allele,
                         term)) %>%
    ggplot(aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
    geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
    geom_label(aes(label = term, fill = estimate),
               family = "Arial", size = 2.2,
               label.padding = unit(0.8, "mm"),
               label.r = unit(0.4, "mm"),
               color = "black", label.size = 0) +
    scale_color_gradient2(high = pastel_red,
                          mid = "grey90",
                          low = pastel_blue) +
    scale_fill_gradient2(high = pastel_red,
                         mid = "grey90",
                         low = pastel_blue) +
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
        "estimated effect of mutation on dependency on *KIAA1257*"
    ))
saveFigRds(kiaa1257_coef_plot, "PAAD_G12R_KIAA1257_coef-plot.rds")
ggsave_wrapper(
    kiaa1257_coef_plot,
    plot_path(GRAPHS_DIR, "PAAD_G12R_KIAA1257_coef-plot.svg"),
    "small"
)


g12r_dnah5_line <- tibble(
    name = c("DNAH5", "reduced comutation", "G12R"),
    color = c("white", comut_updown_pal[["reduced"]], "white"),
    y = c(1, 1.7, 1)
) %>%
    mutate(name = fct_inorder(name)) %>%
    ggplot(aes(x = name, y = y)) +
    geom_segment(x = "DNAH5", xend = "G12R", y = 1, yend = 1,
                 color = comut_updown_pal[["reduced"]],
                 size = 1.4) +
    geom_label(aes(label = name, fill = name, color = color),
               family = "Arial", size = 2.2,
               label.r = unit(0.4, "mm"),
               label.size = 0) +
    scale_color_identity() +
    scale_fill_manual(values = kiaa1257_pal, guide = FALSE) +
    scale_y_continuous(limits = c(0.5, 1.9), expand = c(0, 0)) +
    theme_void(base_size = 7, base_family = "Arial")
saveFigRds(g12r_dnah5_line, "PAAD_G12R_KIAA1257_G12R-DNAH5_line-plot.rds")


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



#### ---- PAAD G12D: FKBP1A ---- ####

# G12D has increased comut with GPR98 and RNF43.
# When one is mutated, dep. on FKBP1A is increased.
# When KRAS and one other are mutated, dep. on FKBP1A is even stronger.

# group: linking

fkbp1a_res <- synlet_comut_model_res %>%
    filter(cancer == "PAAD" & allele == "G12D" & hugo_symbol == "FKBP1A")

fkbp1a_comuts <- c("GPR98", "RNF43")


fkbp1a_coef_plot <- fkbp1a_res$fit[[1]]$elastic_model %>%
    broom::tidy() %>%
    filter(term != "(Intercept)") %>%
    mutate(
        term = str_replace(term, "kras_allele", fkbp1a_res$fit[[1]]$allele),
        term = str_replace(term, ":", " & "),
        term = fct_reorder(term, estimate)
    ) %>%
    ggplot(aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
    geom_segment(aes(yend = term, color = estimate), xend = 0, size = 1) +
    geom_label(aes(label = term, fill = estimate),
               family = "Arial", size = 2.2,
               label.padding = unit(0.8, "mm"),
               label.r = unit(0.4, "mm"),
               color = "black", label.size = 0) +
    scale_color_gradient2(high = pastel_red,
                          mid = "grey90",
                          low = pastel_blue) +
    scale_fill_gradient2(high = pastel_red,
                         mid = "grey90",
                         low = pastel_blue) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.1))) +
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
        "estimated effect of mutation on dependency on *FKBP1A*"
    ))
saveFigRds(fkbp1a_coef_plot, "PAAD_G12D_FKBP1A_coef-plot.rds")

fkbp1a_box_plot <- fkbp1a_res$fit[[1]]$data %>%
    mutate(
        color = case_when(
            GPR98 + RNF43 == 2 ~ "GPR98 & RNF43",
            GPR98 == 1 ~ "GPR98",
            RNF43 == 1 ~ "RNF43",
            kras_allele == 1 ~ "G12D",
            TRUE ~ "WT"),
        shape = ifelse(kras_allele == 1, "G12D", "other"),
        group = ifelse(color == "WT", "WT", "mutant"),
        group = factor(group, levels = c("WT", "mutant"))
    ) %>%
    ggplot(aes(x = group, y = gene_effect)) +
    geom_hline(yintercept = 0, size = 0.2, color = "grey30") +
    ggforce::geom_mark_hull(
        aes(filter = gene_effect < -0.18,
            label = "comutations with strongest effects"),
        fill = "grey75",
        expand = unit(3, "mm"),
        radius = unit(1, "mm"),
        size = 0.3, colour = NA,
        label.family = "Arial",
        label.fontsize = 6,
        label.fontface = "plain",
        label.fill = NA,
        label.colour = "black",
        label.margin = margin(1, 1, 1, 1, "mm"),
        label.buffer = unit(5, "mm"),
        con.colour = "grey50",
        con.size = 0.3
    ) +
    ggbeeswarm::geom_quasirandom(
        aes(color = color, shape = shape),
        size = 2
    ) +
    scale_color_manual(
        values = fkbp1a_box_pal,
        guide = guide_legend(
            direction = "vertical",
            nrow = 2,
            order = 20
        )
    ) +
    scale_shape_manual(
        values = c(17, 16),
        guide = guide_legend(
            direction = "vertical",
            nrow = 2,
            order = 20
        )
    ) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.title = element_markdown(),
        legend.spacing.y = unit(1, "mm"),
        legend.key.height = unit(3, "mm")
    ) +
    labs(y = glue("dependency score of *{fkbp1a_res$hugo_symbol[[1]]}*"),
         color = "mutations",
         shape = "*KRAS* allele")
saveFigRds(fkbp1a_box_plot, "PAAD_G12D_FKBP1A_box-plot.rds")


fkbp1a_genetic_graph <- tribble(
    ~ from,   ~ to,     ~ edge_name,  ~ interaction_type, ~ interaction_spec,
    "G12D",   "GPR98",  "",           "comutation",       "increased comut.",
    "G12D",   "RNF43",  "",           "comutation",       "increased comut.",
    "G12D",   "RYR1",   "",           "comutation",       "reduced comut.",
    "G12D",   "FKBP1A", "",           "dependency",       "more dep.",
    "FKBP1A", "RYR1",   "calcium\nregulation", "PPI",              "PPI regulation",
    "RNF43",  "FKBP1A",  "",          "dependency",       "more dep.",
    "GPR98",  "FKBP1A",  "",          "dependency",       "more dep.",
) %>%
    as_tbl_graph() %E>%
    mutate(interaction_spec = fct_inorder(interaction_spec))



fkbp1a_graph_layout <- ggraph::create_layout(fkbp1a_genetic_graph, "nicely")
fkbp1a_graph_layout_attrs <- attributes(fkbp1a_graph_layout)

fkbp1a_graph_layout$x <- c(1, 2, 3, 3, 3)
fkbp1a_graph_layout$y <- c(2, 4, 2, 1, 3)

fkbp1a_graph_plot <- ggraph(fkbp1a_graph_layout) +
    geom_edge_diagonal(
        aes(edge_color = interaction_spec),
        strength = 0.5,
        width = 0.8
        # start_cap = circle(3.5, "mm"),
        # end_cap = circle(3.5, "mm")
    ) +
    geom_node_label(
        aes(label = name),
        family = "Arial", size = 2.2,
        label.padding = unit(0.8, "mm"),
        label.r = unit(0.4, "mm"),
        color = "black", label.size = 0,
        fill = "grey90"
    ) +
    scale_edge_color_manual(
        values = fkbp1a_edge_pal,
        guide = guide_legend(
            title = "interaction type",
            title.position = "top",
            nrow = 1,
            order = 10
        )
    ) +
    theme_graph(base_size = 7, base_family = "Arial") +
    theme(
        legend.position = "bottom",
        legend.title = element_markdown(face = "bold"),
        legend.spacing.y = unit(1, "mm"),
        legend.key.height = unit(3, "mm")
    )
saveFigRds(fkbp1a_graph_plot, "PAAD_G12D_FKBP1A_genetic-graph-plot.rds")
ggsave_wrapper(
    fkbp1a_graph_plot,
    plot_path(GRAPHS_DIR, "PAAD_G12D_FKBP1A_genetic-graph.svg"),
    "small"
)


fkbp1a_patch <- (
    (fkbp1a_coef_plot / fkbp1a_box_plot) + plot_layout(heights = c(1, 2)) |
    wrap_elements(full = fkbp1a_graph_plot)
)
ggsave_wrapper(
    fkbp1a_patch,
    plot_path(GRAPHS_DIR, "PAAD_G12D_FKBP1A_patch.svg"),
    "wide"
)
