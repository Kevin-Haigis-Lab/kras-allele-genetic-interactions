# Visualizing the genetic dependency connections.

GRAPHS_DIR <- "10_18_alluvial-plots-comparing-dependency"
reset_graph_directory(GRAPHS_DIR)

library(ggalluvial)


prepare_depdency_data_for_alluvial <- function(data) {
    mod_data <- data %>%
        mutate(
            gene_cls = factor(gene_cls),
            avg_gene_effect_rnd = cut(avg_gene_effect,
                                      breaks = seq(-4, 4, 0.5)),
            avg_gene_effect_rnd = fct_rev(avg_gene_effect_rnd)
        )
    return(mod_data)
}


dependency_cancer_alluvial_plot <- function(data) {
    p <- data %>%
        ggplot(aes(x = kras_allele,
                   alluvium = hugo_symbol,
                   stratum = avg_gene_effect_rnd)) +
        ggalluvial::geom_alluvium(aes(fill = gene_cls)) +
        ggalluvial::geom_stratum() +
        geom_text(stat = "stratum", infer.label = TRUE, size = 2)
    invisible(p)
}


style_alluvial_plot <- function(p, cancer) {
    p <- p +
        scale_fill_viridis_d() +
        scale_x_discrete(expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw(base_size = 7, base_family = "arial") +
        theme(
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        ) +
        labs(
            x = "KRAS allele",
            y = "average gene effect",
            title = str_to_upper(cancer)
        )
    invisible(p)
}


save_alluvial_plot <- function(p, cancer, template) {
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, glue(template)),
        "wide"
    )
    return(NULL)
}


save_alluvial_patch <- function(ps) {
    patch <- wrap_plots(ps, ncol = 1)
    ggsave_wrapper(
        patch,
        plot_path(GRAPHS_DIR, "genedep-alluvial_patchwork.svg"),
        "tall"
    )
    return(NULL)
}


depmap_model_workflow_res %>%
    filter_depmap_model_workflow_res() %>%
    select(cancer, hugo_symbol, data) %>%
    unnest(data) %>%
    group_by(cancer, hugo_symbol) %>%
    mutate(gene_effect = scale(gene_effect)[, 1]) %>%
    group_by(cancer, hugo_symbol, kras_allele) %>%
    summarise(avg_gene_effect = mean(gene_effect)) %>%
    ungroup() %>%
    left_join(depmap_gene_clusters, by = c("cancer", "hugo_symbol")) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        data = purrr::map(data, prepare_depdency_data_for_alluvial),
        alluvial_plot = purrr::map(data, dependency_cancer_alluvial_plot),
        alluvial_plot = map2(alluvial_plot, cancer, style_alluvial_plot),
        save_plot = map2(alluvial_plot, cancer, save_alluvial_plot,
                         template = "gendep-alluvial_{cancer}.svg"),
        save_patch = save_alluvial_patch(alluvial_plot)
    )
