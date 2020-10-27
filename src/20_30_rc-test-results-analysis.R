# Analyze the results of the RC-test for comutation and mutual exclusivity

library(tidygraph)
library(ggraph)

#### ---- Barplots for number of sig. results ---- ####

info(logger, "Making barplots of general results of RC-test.")

sig_results_barplot <- rc_test_results %>%
  filter(num_mut_per_cancer >= 3) %>%
  filter(p_val < 0.05 & t_AM >= 3) %>%
  ggplot(aes(x = allele)) +
  facet_wrap(~cancer, scales = "free") +
  geom_bar(aes(fill = cancer)) +
  scale_fill_manual(values = cancer_palette, guide = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic() +
  theme(
    text = element_text(family = "arial"),
    axis.text.x = element_text(angle = 30, hjust = 1.0)
  )


save_path <- plot_path("20_30_rc-test-results-analysis", "sig_results_barplot.svg")
info(logger, glue("Making plot {save_path}."))
ggsave_wrapper(
  sig_results_barplot,
  save_path,
  "medium"
)


sig_results_barplot_septest <- rc_test_results %>%
  filter(num_mut_per_cancer >= 3) %>%
  filter(p_val < 0.05 & t_AM >= 3) %>%
  ggplot(aes(x = allele)) +
  facet_wrap(~cancer, scales = "free") +
  geom_bar(aes(fill = rc_test_type), position = "dodge") +
  scale_fill_manual(
    values = comut_mutex_pal
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic() +
  theme(
    text = element_text(family = "arial"),
    axis.text.x = element_text(angle = 30, hjust = 1.0),
    legend.position = "bottom"
  )

save_path <- plot_path(
  "20_30_rc-test-results-analysis",
  "sig_results_barplot_septest.svg"
)
info(logger, glue("Making plot {save_path}."))
ggsave_wrapper(sig_results_barplot_septest, save_path, "medium")



#### ---- Networks ---- ####

rc_test_gr <- rc_test_results %>%
  filter(num_mut_per_cancer >= 10) %>%
  filter(p_val < 0.01 & t_AM >= 3) %>%
  mutate(allele = paste("KRAS", allele)) %>%
  select(allele, hugo_symbol, cancer, rc_test_type, p_val, t_AM) %>%
  as_tbl_graph() %N>%
  mutate(is_kras = str_detect(name, "KRAS "))


for (CANCER in unique(rc_test_results$cancer)) {
  rc_gr_plot <- rc_test_gr %E>%
    filter(cancer == !!CANCER) %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(color_label = ifelse(is_kras, name, NA)) %>%
    ggraph(layout = "stress") +
    geom_edge_link(
      aes(
        color = rc_test_type,
        width = -log(p_val)
      ),
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

  save_path <- plot_path(
    "20_30_rc-test-results-analysis",
    glue("{CANCER}_comut-mutex_network.svg")
  )
  info(logger, glue("Making plot {save_path}."))
  ggsave_wrapper(rc_gr_plot, save_path, "large")
}
