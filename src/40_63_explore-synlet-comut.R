# Explore some of the results from the modeling of synthetic lethality
# using comutation interactions.

GRAPHS_DIR <- "40_63_explore-synlet-comut"
reset_graph_directory(GRAPHS_DIR)

pastel_red <- "#FF8FAC"
pastel_blue <- "#4096B3"

set.seed(0)


#### ---- MASKING ---- ####

# Make the coefficient plot for the masking group.
make_mask_coef_plot <- function(cancer, allele, hugo_symbol, fit, ...) {

    title <- paste0(
        cancer, "<br>",
        glue("estimated effect of mutation on dependency on *{hugo_symbol}*")
    )

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
                   label.padding = unit(0.4, "mm"),
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
            plot.title = element_markdown(hjust = 0),
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(title = title)
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
        geom_boxplot(fill = NA, outlier.shape = NA) +
        ggbeeswarm::geom_quasirandom(size = 1, alpha = 0.9, groupOnX = TRUE) +
        scale_color_manual(values = mask_pal, drop = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            legend.position = "none",
            axis.title.y = element_markdown(),
            axis.title.x = element_blank()
        ) +
        labs(y = glue("dependency score of *{hugo_symbol}*"))

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


# A standard coefficient plot for the additive group of models.
make_additive_coef_plot <- function(cancer, allele, hugo_symbol, fit, coeffs,
                                    ...) {

    title <- paste0(
        cancer, "<br>",
        "estimated effect of mutation<br>",
        glue("on dependency on *{hugo_symbol}*")
    )

    p <- broom::tidy(fit$elastic_model) %>%
        filter(term != "(Intercept)") %>%
        filter(term %in% !!coeffs) %>%
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
            plot.title = element_markdown(hjust = 0),
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(title = title)
    plt_name <- glue("{cancer}_{allele}_{hugo_symbol}_coef-plot.svg") %>%
        as.character()
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, plt_name),
        "small"
    )
    return(p)
}


# A collection of standard value for the plots for additive models.
stds <- list(
    dotplot = list(size = 3, alpha = 0.8),
    line = list(grey = "grey80", alpha = 1, size = 0.7,
                linejoin = "round", lineend = "round", lty = 1),
    linept = list(grey = "grey80", size = 0.7, alpha = 1),
    zerohline = list(size = 0.2, grey = "grey30"),
    gr_edge = list(width = 1, label_size = 2.2, vjust = -0.5,
                   lineheight = 0.8),
    gr_node = list(size = 2.2, label.padding = unit(0.8, "mm"),
                   label.r = unit(0.4, "mm"), label.size = 0,
                   grey = "grey85"),
    pal = list(neither = "grey70", other = "grey30")
)



#### ---- COAD G12D: SRSF5 ---- ####

# G12D has reduced comut with HECW1 and increased comut with APC.
# When all three are mutated, there is less dep. on SRSF5.

# group: linking

srsf5_res <- synlet_comut_model_res %>%
    filter(cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5")

srsf5_comuts <- c("APC", "HECW1")
srsf5_terms <- c("kras_allele", srsf5_comuts, "kras_allele:APC")

plots <- srsf5_res %>%
    pmap(make_additive_coef_plot, coeffs = srsf5_terms)
srsf5_coef_plot <- plots[[1]]
saveFigRds(srsf5_coef_plot, "COAD_G12D_SRSF5_coef-plot.rds")


srsf_box_pal <- c(short_allele_pal["G12D"],
                  "APC/HECW1" = stds$pal$other,
                  neither = stds$pal$neither)

srsf5_box_plot_data <- srsf5_res$fit[[1]]$data %>%
    mutate(
        color = case_when(
            kras_allele == 1 ~ !!srsf5_res$allele,
            APC + HECW1 > 0 ~ "APC/HECW1",
            TRUE ~ "neither"),
        color = factor(color,
                       levels = c(srsf5_res$allele, "APC/HECW1", "neither")),
        mutation = case_when(
            APC + HECW1 == 2 ~ "APC & HECW1",
            APC == 1 ~ "APC",
            HECW1 == 1 ~ "HECW1",
            TRUE ~ "neither"),
        mutation = factor(mutation,
                          levels = c("neither", "APC", "HECW1", "APC & HECW1"))
    )

srsf5_box_plot_summary <- srsf5_box_plot_data %>%
    group_by(mutation) %>%
    summarise(gene_effect = mean(gene_effect)) %>%
    ungroup()

srsf5_box_plot <- srsf5_box_plot_data %>%
    ggplot(aes(mutation, gene_effect)) +
    geom_hline(yintercept = 0,
               size = stds$zerohline$size,
               color = stds$zerohline$grey) +
    geom_line(aes(group = "a"),
              data = srsf5_box_plot_summary,
              size = stds$line$size,
              color = stds$line$grey,
              lty = stds$line$lty,
              linejoin = stds$line$linejoin,
              lineend = stds$line$lineend) +
    geom_point(data = srsf5_box_plot_summary,
               size = stds$linept$size,
               color = stds$linept$grey) +
    ggbeeswarm::geom_quasirandom(aes(color = color),
                                 size = stds$dotplot$size,
                                 alpha = stds$dotplot$alpha) +
    scale_color_manual(values = srsf_box_pal, drop = TRUE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.ticks = element_blank(),
        legend.title = element_markdown()
    ) +
    labs(y = glue("dependency score of *{srsf5_res$hugo_symbol[[1]]}*"),
         color = "mutation")
saveFigRds(srsf5_box_plot, "COAD_G12D_SRSF5_box-plot.rds")


srsf5_comut_gr <- genetic_interaction_gr %E>%
    filter(cancer == "COAD") %>%
    mutate(genetic_interaction = switch_comut_terms(genetic_interaction),
           edge_label = paste0(genetic_interaction, "\ncomutation")) %N>%
    mutate(name = str_remove(name, "KRAS_")) %>%
    filter(name %in% c("G12D", srsf5_comuts)) %>%
    mutate(node_text_color = ifelse(name %in% kras_dark_lbls, "white", "black"))

srsf5_comut_gr_layout <- create_layout(srsf5_comut_gr, layout = "nicely")
srsf5_comut_gr_layout$x <- c(2, 3, 1)
srsf5_comut_gr_layout$y <- rep(0, 3)
srsf5_comut_plot <- ggraph(srsf5_comut_gr_layout) +
    geom_edge_link(aes(color = genetic_interaction,
                       label = edge_label),
                   vjust = stds$gr_edge$vjust,
                   label_size = stds$gr_edge$label_size,
                   family = "Arial",
                   width = stds$gr_edge$width,
                   lineheight = stds$gr_edge$lineheight) +
    geom_node_label(aes(label = name, fill = name, color = node_text_color),
                    family = "Arial",
                    size = stds$gr_node$size,
                    label.padding = stds$gr_node$label.padding,
                    label.r = stds$gr_node$label.r,
                    label.size = stds$gr_node$label.size) +
    scale_edge_color_manual(values = comut_updown_pal, guide = FALSE) +
    scale_color_identity(guide = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    scale_fill_manual(values = short_allele_pal,
                      na.value = stds$gr_node$grey,
                      guide = FALSE) +
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

kiaa1257_box_plot_data <- kiaa1257_res$fit[[1]]$data %>%
    mutate(
        mutation = case_when(
            kras_allele + DNAH5 == 2 ~ "DNAH5 & G12R",
            kras_allele == 1 ~ "G12R",
            DNAH5 == 1 ~ "DNAH5",
            TRUE ~ "neither"
        ),
        mutation = factor(
            mutation,
            levels = c("neither", "DNAH5", "G12R", "DNAH5 & G12R")
        )
    )

kiaa1257_box_plot_summary <- kiaa1257_box_plot_data %>%
    group_by(mutation) %>%
    summarise(gene_effect = mean(gene_effect)) %>%
    ungroup()

kiaa1257_box_pal <- c(
    "neither" = "grey70",
    "DNAH5" = "grey30",
    "G12R" = short_allele_pal[["G12R"]],
    "DNAH5 & G12R" = "#5b2280"
)

kiaa1257_box_plot <- kiaa1257_box_plot_data %>%
    ggplot(aes(mutation, gene_effect)) +
    geom_hline(yintercept = 0,
               size = stds$zerohline$size,
               color = stds$zerohline$grey) +
    geom_line(aes(group = "a"),
              data = kiaa1257_box_plot_summary,
              size = stds$line$size,
              color = stds$line$grey,
              lty = stds$line$lty,
              linejoin = stds$line$linejoin,
              lineend = stds$line$lineend) +
    geom_point(data = kiaa1257_box_plot_summary,
               size = stds$linept$size,
               color = stds$linept$grey) +
    ggbeeswarm::geom_quasirandom(aes(color = mutation),
                                 size = stds$dotplot$size,
                                 alpha = stds$dotplot$alpha) +
    scale_color_manual(values = kiaa1257_box_pal, drop = TRUE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.title = element_blank(),
        legend.position = "right",
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


plots <- kiaa1257_res %>%
    pmap(make_additive_coef_plot, coeffs = c("DNAH5", "kras_allele"))
kiaa1257_coef_plot <- plots[[1]]
saveFigRds(kiaa1257_coef_plot, "PAAD_G12R_KIAA1257_coef-plot.rds")
ggsave_wrapper(
    kiaa1257_coef_plot,
    plot_path(GRAPHS_DIR, "PAAD_G12R_KIAA1257_coef-plot.svg"),
    "small"
)


kiaa1257_comut_gr <- genetic_interaction_gr %E>%
    filter(cancer == "PAAD") %>%
    mutate(genetic_interaction = switch_comut_terms(genetic_interaction),
           edge_label = paste0(genetic_interaction, "\ncomutation")) %N>%
    mutate(name = str_remove(name, "KRAS_")) %>%
    filter(name %in% c("G12R", "DNAH5")) %>%
    mutate(node_text_color = ifelse(name %in% kras_dark_lbls, "white", "black"))

kiaa1257_comut_gr_layout <- create_layout(kiaa1257_comut_gr, layout = "nicely")
kiaa1257_comut_gr_layout$x <- c(2, 1)
kiaa1257_comut_gr_layout$y <- rep(0, 2)
kiaa1257_comut_plot <- ggraph(kiaa1257_comut_gr_layout) +
    geom_edge_link(aes(color = genetic_interaction,
                       label = edge_label),
                   vjust = stds$gr_edge$vjust,
                   label_size = stds$gr_edge$label_size,
                   family = "Arial",
                   width = stds$gr_edge$width,
                   lineheight = stds$gr_edge$lineheight) +
    geom_node_label(aes(label = name, fill = name, color = node_text_color),
                    family = "Arial",
                    size = stds$gr_node$size,
                    label.padding = stds$gr_node$label.padding,
                    label.r = stds$gr_node$label.r,
                    label.size = stds$gr_node$label.size) +
    scale_edge_color_manual(values = comut_updown_pal, guide = FALSE) +
    scale_color_identity(guide = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.4, 0.4))) +
    scale_fill_manual(values = short_allele_pal,
                      na.value = stds$gr_node$grey,
                      guide = FALSE) +    theme_graph() +
    theme(
        plot.margin = margin(3, 0, 0, 0, unit = "mm")
    )
saveFigRds(kiaa1257_comut_plot, "PAAD_G12R_KIAA1257_comut-graph-plot.rds")


kiaa1257_patch <- (kiaa1257_coef_plot / kiaa1257_box_plot / kiaa1257_comut_plot) +
    plot_layout(heights = c(3, 9, 1))
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

plots <- fkbp1a_res %>%
    pmap(make_additive_coef_plot,
         coeffs = c("RNF43", "kras_allele:GPR98", "kras_allele"))
fkbp1a_coef_plot <- plots[[1]]
saveFigRds(fkbp1a_coef_plot, "PAAD_G12D_FKBP1A_coef-plot.rds")

fkbp1a_box_pal <- c(
    "neither" = "grey70",
    "G12D" = short_allele_pal[["G12D"]],
    "GPR98/RNF43" = "grey30"
)

fkbp1a_box_plot_data <- fkbp1a_res$fit[[1]]$data %>%
    mutate(
        mutation = case_when(
            GPR98 + RNF43 == 2 ~ "GPR98 & RNF43",
            GPR98 == 1 ~ "GPR98",
            RNF43 == 1 ~ "RNF43",
            TRUE ~ "neither"
        ),
        color = case_when(
            kras_allele == 1 ~ "G12D",
            GPR98 + RNF43 >= 1 ~ "GPR98/RNF43",
            TRUE ~ "neither"
        ),
        mutation = factor(
            mutation,
            levels = c("neither", "GPR98", "RNF43", "GPR98 & RNF43")
        ),
        color = factor(
            color,
            levels = c("neither", "GPR98/RNF43", "G12D")
        )
    )

fkbp1a_box_plot_summary <- fkbp1a_box_plot_data %>%
    group_by(kras_allele, mutation) %>%
    summarise(gene_effect = mean(gene_effect)) %>%
    ungroup() %>%
    arrange(mutation)


fkbp1a_box_plot <- fkbp1a_box_plot_data %>%
    ggplot(aes(x = mutation, y = gene_effect)) +
    geom_hline(yintercept = 0,
               size = stds$zerohline$size,
               color = stds$zerohline$grey) +
    geom_line(data = filter(fkbp1a_box_plot_summary, kras_allele == 0),
              group = "a",
              color = stds$line$grey,
              size = stds$line$size) +
    geom_point(data = filter(fkbp1a_box_plot_summary, kras_allele == 0),
               group = "a",
               color = stds$linept$grey,
               size = stds$linept$size) +
    geom_line(data = filter(fkbp1a_box_plot_summary, kras_allele == 1),
              group = "a",
              color = "#87a7ff",
              size = stds$line$size) +
    geom_point(data = filter(fkbp1a_box_plot_summary, kras_allele == 1),
               group = "a",
               color = "#87a7ff",
               size = stds$linept$size) +
    geom_point(aes(color = color),
               position = position_jitterdodge(dodge.width = 0.4),
               size = stds$dotplot$size,
               alpha = stds$dotplot$alpha) +
    scale_color_manual(values = fkbp1a_box_pal, drop = TRUE) +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_blank()
    ) +
    labs(y = glue("dependency score of *{fkbp1a_res$hugo_symbol[[1]]}*"),
         color = "*KRAS* allele")
saveFigRds(fkbp1a_box_plot, "PAAD_G12D_FKBP1A_box-plot.rds")


fkbp1a_comut_gr <- genetic_interaction_gr %E>%
    filter(cancer == "PAAD") %>%
    mutate(genetic_interaction = switch_comut_terms(genetic_interaction),
           edge_label = paste0(genetic_interaction, "\ncomutation")) %N>%
    mutate(name = str_remove(name, "KRAS_")) %>%
    filter(name %in% c("G12D", fkbp1a_comuts)) %>%
    mutate(node_text_color = ifelse(name %in% kras_dark_lbls, "white", "black"))

fkbp1a_comut_gr_layout <- create_layout(fkbp1a_comut_gr, layout = "nicely")
fkbp1a_comut_gr_layout$x <- c(2, 3, 1)
fkbp1a_comut_gr_layout$y <- rep(0, 3)
fkbp1a_comut_plot <- ggraph(fkbp1a_comut_gr_layout) +
    geom_edge_link(aes(color = genetic_interaction,
                       label = edge_label),
                   vjust = stds$gr_edge$vjust,
                   label_size = stds$gr_edge$label_size,
                   family = "Arial",
                   width = stds$gr_edge$width,
                   lineheight = stds$gr_edge$lineheight) +
    geom_node_label(aes(label = name, fill = name, color = node_text_color),
                    family = "Arial",
                    size = stds$gr_node$size,
                    label.padding = stds$gr_node$label.padding,
                    label.r = stds$gr_node$label.r,
                    label.size = stds$gr_node$label.size) +
    scale_edge_color_manual(values = comut_updown_pal, guide = FALSE) +
    scale_color_identity(guide = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    scale_fill_manual(values = short_allele_pal,
                      na.value = stds$gr_node$grey,
                      guide = FALSE) +
    theme_graph() +
    theme(
        plot.margin = margin(3, 0, 0, 0, unit = "mm")
    )
saveFigRds(fkbp1a_comut_plot, "PAAD_G12D_FKBP1A_comut-graph-plot.rds")


fkbp1a_patch <- (fkbp1a_coef_plot / fkbp1a_box_plot / fkbp1a_comut_plot) +
    plot_layout(heights = c(3, 9, 1))
ggsave_wrapper(
    fkbp1a_patch,
    plot_path(GRAPHS_DIR, "PAAD_G12D_FKBP1A-patch.svg"),
    "small"
)



#### ---- Comutation graphs for masking relationships ---- ####


masking_comut_graph <- function(cancer, allele, gene) {
    comut_gr <- genetic_interaction_gr %E>%
        filter(cancer == !!cancer) %>%
        mutate(genetic_interaction = switch_comut_terms(genetic_interaction),
               edge_label = paste0(genetic_interaction, "\ncomutation")) %N>%
        mutate(name = str_remove(name, "KRAS_")) %>%
        filter(name %in% c(allele, gene)) %>%
        mutate(node_text_color = ifelse(name %in% kras_dark_lbls,
                                        "white", "black"))

    comut_gr_layout <- create_layout(comut_gr, layout = "nicely")
    comut_gr_layout$x <- c(2, 1)
    comut_gr_layout$y <- rep(0, 2)
    comut_plot <- ggraph(comut_gr_layout) +
        geom_edge_link(aes(color = genetic_interaction,
                           label = edge_label),
                       vjust = stds$gr_edge$vjust,
                       label_size = stds$gr_edge$label_size,
                       family = "Arial",
                       width = stds$gr_edge$width,
                       lineheight = stds$gr_edge$lineheight) +
        geom_node_label(aes(label = name, fill = name, color = node_text_color),
                        family = "Arial",
                        size = stds$gr_node$size,
                        label.padding = stds$gr_node$label.padding,
                        label.r = stds$gr_node$label.r,
                        label.size = stds$gr_node$label.size) +
        scale_edge_color_manual(values = comut_updown_pal, guide = FALSE) +
        scale_color_identity(guide = FALSE) +
        scale_x_continuous(expand = expansion(mult = c(0.4, 0.4))) +
        scale_fill_manual(values = short_allele_pal,
                          na.value = stds$gr_node$grey,
                          guide = FALSE) +    theme_graph() +
        theme(
            plot.margin = margin(3, 0, 0, 0, unit = "mm")
        )
    return(comut_plot)
}

coad_tp53_comut_graph <- masking_comut_graph("COAD", "G12D", "TP53")
paad_smad4_comut_graph <- masking_comut_graph("PAAD", "G12D", "SMAD4")

saveFigRds(coad_tp53_comut_graph, "COAD_G12D_TP53_comut-graph-plot.rds")
saveFigRds(paad_smad4_comut_graph, "PAAD_G12D_SMAD4_comut-graph-plot.rds")

