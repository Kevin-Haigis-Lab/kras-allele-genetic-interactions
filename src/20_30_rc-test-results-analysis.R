# Analyze the results of the RC-test for comutation and mutual exclusivity

library(tidygraph)
library(ggraph)

#### ---- Barplots for number of sig. results ---- ####


sig_results_barplot <- rc_test_results %>%
    filter(p_val < 0.05 & t_AM > 2) %>%
    ggplot(aes(x = allele)) +
    facet_wrap(~cancer, scales = "free") +
    geom_bar(aes(fill = cancer)) +
    scale_fill_manual(values = cancer_palette, guide = FALSE) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    theme_classic() +
    theme(
        text = element_text(family = "arial"),
        axis.text.x = element_text(angle = 30, hjust = 1.0)
    )

ggsave_wrapper(
    sig_results_barplot,
    plot_path("20_30_rc-test-results-analysis", "sig_results_barplot.svg"),
    "medium")


sig_results_barplot_septest <- rc_test_results %>%
    filter(p_val < 0.05 & t_AM > 2) %>%
    ggplot(aes(x = allele)) +
    facet_wrap(~cancer, scales = "free") +
    geom_bar(aes(fill = test), position = "dodge") +
    scale_fill_manual(
        values = comut_mutex_pal) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    theme_classic() +
    theme(
        text = element_text(family = "arial"),
        axis.text.x = element_text(angle = 30, hjust = 1.0),
        legend.position = "bottom"
    )

ggsave_wrapper(
    sig_results_barplot_septest,
    plot_path("20_30_rc-test-results-analysis", "sig_results_barplot_septest.svg"),
    "medium")


#### ---- Networks ---- ####

rc_test_gr <- rc_test_results %>%
    filter(p_val < 0.05 & t_AM > 2) %>%
    mutate(allele = paste("KRAS", allele)) %>%
    select(allele, hugo_symbol, cancer, test, p_val, t_AM) %>%
    as_tbl_graph() %N>%
    mutate(is_kras = str_detect(name, "KRAS "))


for (CANCER in unique(rc_test_results$cancer)) {
    rc_gr_plot <- rc_test_gr %E>%
        filter(cancer == !!CANCER) %N>%
        filter(centrality_degree(mode = "all") > 0) %>%
        mutate(color_label = ifelse(is_kras, name, NA)) %>%
        ggraph(layout = "stress") +
        geom_edge_link(
            aes(color = test,
                width = -log(p_val + 0.000001)),
            alpha = 0.7
        ) +
        scale_edge_width_continuous(range = c(0.1, 1)) +
        scale_edge_color_manual(values = comut_mutex_pal) +
        geom_node_point(
            aes(color = color_label),
            size = 1
        ) +
        geom_node_text(
            aes(label = color_label),
            repel = TRUE,
            family = "arial"
        ) +
        scale_color_manual(values = allele_palette, na.value = NA) +
        theme_graph() +
        theme(
            text = element_text(family = "arial")
        ) +
        labs(
            title = glue("{CANCER} KRAS allele-specific genetic interaction network"),
            edge_color = "interaction",
            edge_width = "-log(p-val.)",
            color = "KRAS allele"
        )
    ggsave_wrapper(
        rc_gr_plot,
        plot_path("20_30_rc-test-results-analysis",
                  glue("{CANCER}_comut-mutex_network.svg")),
        "large")
}
