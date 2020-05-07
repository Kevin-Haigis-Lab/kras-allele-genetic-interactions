# PMID: 26353932, 21169473

# KRAS G12C linked to reduced response to cisplatin treatment due to faster or
# better BER. This aligns with our findings on the increased dependency on
# the Fanonci anemia pathway (resolves interstrand crosslinks).

GRAPHS_DIR <- "10_56_luad_depmap_G12C-cisplatin-sensitivity"
reset_graph_directory(GRAPHS_DIR)


# depmap_model_workflow_res %>%
#     filter(cancer == "LUAD" & hugo_symbol == "POLB") %>%
#     pull(allele_aov) %>%
#     purrr::map(~ broom::glance(.x)) %>%
#     knitr::kable()

# polb_boxplot <- model1_tib %>%
#     filter(cancer == "LUAD" & hugo_symbol == "POLB") %>%
#     select(cancer, hugo_symbol, data) %>%
#     unnest(data) %>%
#     ggplot(aes(x = allele, y = gene_effect)) +
#     geom_boxplot(aes(color = allele), outlier.shape = NA) +
#     geom_jitter(aes(color = allele), width = 0.2) +
#     scale_color_manual(values = short_allele_pal) +
#     theme_bw(base_size = 7, base_family = "arial") +
#     theme(
#         axis.title.x = element_blank()
#     ) +
#     labs(
#         y = "dependency score"
#     )
# ggsave_wrapper(
#     polb_boxplot,
#     plot_path(GRAPHS_DIR, "polb_boxplot.svg"),
#     "small"
# )


# genetic_interaction_df %>%
#     filter(cancer == "LUAD" & hugo_symbol == "POLB")
