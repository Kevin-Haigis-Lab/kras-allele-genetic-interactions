# Interesting comutation of STK11 with G12C in LUAD.
# Mutation D194Y is more frequent.

GRAPHS_DIR <- "20_70_luad-g12c-stk11"
reset_graph_directory(GRAPHS_DIR)


cancer_coding_av_muts_df %>%
    filter(cancer == "LUAD" & ras_allele == "KRAS_G12C") %>%
    filter(hugo_symbol == "STK11") %>%
    filter(mutation_type == "missense_mutation") %>%
    filter(amino_acid_change == "D194Y") %>%
    select(sift_pred, polyphen2_hdiv_pred, polyphen2_hvar_pred, lrt_pred,
           mutation_taster_pred, mutation_assessor_pred, fathmm_pred,
           meta_svm_pred, meta_lr_pred, clinsig) %>%
    unique() %>%
    unlist()

ccle_mutations %>%
    filter(hugo_symbol == "STK11" & protein_change == "p.D194Y")

model1_tib %>%
    filter(hugo_symbol == "STK11" & cancer == "LUAD") %>%
    select(cancer, hugo_symbol, data) %>%
    unnest(data) %>%
    filter(dep_map_id == "ACH-000698")



stk11_mutations <- cancer_coding_av_muts_df %>%
    filter(cancer == "LUAD" & hugo_symbol == "STK11") %>%
    mutate(group = ifelse(ras_allele == "KRAS_G12C", "G12C", "rest")) %>%
    group_by(group) %>%
    mutate(
        num_samples = n_distinct(tumor_sample_barcode),
        num_stk11_mutations = n()
    ) %>%
    group_by(group, amino_position, mutation_type,
             num_samples, num_stk11_mutations) %>%
    summarise(
        num_mutation_type = n(),
        percent_of_mutations = num_mutation_type / unique(num_stk11_mutations)
    ) %>%
    ungroup()

stk11_mutations %<>%
    mutate(mutation_type = ifelse(mutation_type == "frameshift_deletion",
                                  "frame_shift_del", mutation_type))

stk11_lollipop <- stk11_mutations %>%
    mutate(
        amino_position = as.numeric(amino_position),
        percent_of_mutations = ifelse(group == "G12C", percent_of_mutations, -percent_of_mutations)
    ) %>%
    ggplot(aes(x = amino_position, y = percent_of_mutations)) +
    geom_col(aes(fill = mutation_type), position = "stack") +
    theme_bw(base_size = 7, base_family = "Arial")
ggsave_wrapper(
    stk11_lollipop,
    plot_path(GRAPHS_DIR, "stk11_lollipop.svg"),
    "wide"
)




#### ---- STK11 Lollipop ---- ####

stk11_pal <- c(
    "N-term" = "goldenrod1",
    "NLS" = "mediumorchid1",
    "Kinase" = "darkturquoise",
    "C-term" = "darksalmon"
)

stk11_scheme <- tibble(
    x = c(1, 38, 38, 43, 43, 49, 49, 309, 309, 433),
    domain = c("N-term", "N-term", "NLS", "NLS", "N-term", "N-term",
               "Kinase", "Kinase", "C-term", "C-term")
)

stk11_labels <- tibble(
    x = c(mean(c(1, 48)), mean(c(49, 309)), mean(c(310, 433))),
    domain = c("N-term", "Kinase", "C-term")
)

stk11_scheme_plot <- stk11_scheme %>%
    ggplot(aes(x = x, y = 0)) +
    geom_line(aes(color = domain), size = 4) +
    geom_text(
        data = stk11_labels,
        aes(x = x, y = 0, label = domain),
        size = 2, family = "Arial", fontface = "bold"
    ) +
    scale_x_continuous(expand = expand_scale(mult = c(0.01, 0.01))) +
    scale_color_manual(values = stk11_pal, guide = FALSE) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    theme(
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(-2, 0, -2, 0, "mm")
    ) +
    labs(
        color = "STK11 structure"
    )
ggsave_wrapper(
    stk11_scheme_plot,
    plot_path(GRAPHS_DIR, "stk11_scheme_plot.svg"),
    "small"
)


stk11_g12c_lollipop <- stk11_mutations %>%
    filter(group == "G12C") %>%
    mutate(amino_position = as.numeric(amino_position),
           mutation_type = format_var_names(mutation_type)) %>%
    ggplot(aes(x = amino_position, y = num_mutation_type)) +
    geom_col(aes(fill = mutation_type), position = "stack") +
    geom_hline(yintercept = 0, size = 0.8, color = "black") +
    scale_fill_manual(values = mod_variant_pal) +
    scale_x_continuous(
        expand = expand_scale(mult = c(0.01, 0.01)),
        breaks = seq(25, 433, 25)
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, -2, 0, "mm"),
        axis.ticks = element_blank(),
        panel.border = element_blank()
    ) +
    labs(
        fill = "mutation type",
        y = "num. mut.\nin G12C samples"
    )
ggsave_wrapper(
    stk11_g12c_lollipop,
    plot_path(GRAPHS_DIR, "stk11_g12c_lollipop.svg"),
    "small"
)

stk11_rest_lollipop <- stk11_mutations %>%
    filter(group != "G12C") %>%
    mutate(amino_position = as.numeric(amino_position),
           mutation_type = format_var_names(mutation_type),
           num_mutation_type = -1 * num_mutation_type) %>%
    ggplot(aes(x = amino_position, y = num_mutation_type)) +
    geom_col(aes(fill = mutation_type), position = "stack") +
    geom_hline(yintercept = 0, size = 0.8, color = "black") +
    scale_fill_manual(values = mod_variant_pal) +
    scale_x_continuous(
        expand = expand_scale(mult = c(0.01, 0.01)),
        breaks = seq(25, 433, 25)
    ) +
    scale_y_continuous(
        expand = expand_scale(mult = c(0.02, 0)),
        labels = abs
    ) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
        plot.margin = margin(-2, 0, 0, 0, "mm"),
        axis.ticks = element_blank(),
        panel.border = element_blank()
    ) +
    labs(
        fill = "mutation type",
        y = "num. mut.\nin other samples",
        x = "amino acid position on STK11"
    )
ggsave_wrapper(
    stk11_rest_lollipop,
    plot_path(GRAPHS_DIR, "stk11_rest_lollipop.svg"),
    "small"
)


stk11_lollipop_patch <- (stk11_g12c_lollipop / stk11_scheme_plot / stk11_rest_lollipop) +
    plot_layout(heights = c(2, 1, 2), guides = "collect") &
    theme(
        legend.position = "bottom"
    )
ggsave_wrapper(
    stk11_lollipop_patch,
    plot_path(GRAPHS_DIR, "stk11_lollipop_patch.svg"),
    "wide"
)

saveRDS(
    stk11_lollipop_patch,
    get_fig_proto_path("stk11_lollipop_patch.rds", 5)
)
