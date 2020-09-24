# Interesting comutation of STK11 with G12C in LUAD.
# Mutation D194Y is more frequent.

GRAPHS_DIR <- "20_70_luad-g12c-stk11"
reset_graph_directory(GRAPHS_DIR)


#### ---- EDA ---- ####

cancer_coding_av_muts_df %>%
  filter(cancer == "LUAD" & ras_allele == "KRAS_G12C") %>%
  filter(hugo_symbol == "STK11") %>%
  filter(mutation_type == "missense_mutation") %>%
  filter(amino_acid_change == "D194Y") %>%
  select(
    sift_pred, polyphen2_hdiv_pred, polyphen2_hvar_pred, lrt_pred,
    mutation_taster_pred, mutation_assessor_pred, fathmm_pred,
    meta_svm_pred, meta_lr_pred, clinsig
  ) %>%
  unique() %>%
  unlist()

ccle_mutations %>%
  filter(hugo_symbol == "STK11" & protein_change == "p.D194Y")

model1_tib %>%
  filter(hugo_symbol == "STK11" & cancer == "LUAD") %>%
  select(cancer, hugo_symbol, data) %>%
  unnest(data) %>%
  filter(dep_map_id == "ACH-000698")


#### ---- Prepare STK11 mutation data ---- ####

parse_sift <- function(x) {
  case_when(
    str_detect(x, "D") ~ TRUE,
    TRUE ~ FALSE
  )
}

parse_polyphen <- function(x) {
  case_when(
    str_detect(x, "D") ~ TRUE,
    str_detect(x, "P") ~ TRUE,
    TRUE ~ FALSE
  )
}

parse_fathmm <- function(x) {
  case_when(
    str_detect(x, "D") ~ TRUE,
    TRUE ~ FALSE
  )
}

parse_clinsig <- function(x) {
  case_when(
    str_detect(x, "Path") ~ TRUE,
    TRUE ~ FALSE
  )
}


stk11_mutations <- cancer_coding_av_muts_df %>%
  filter(cancer == "LUAD") %>%
  mutate(
    group = ifelse(ras_allele == "KRAS_G12C", "G12C", "rest"),
    num_luad_samples = n_distinct(tumor_sample_barcode)
  ) %>%
  filter(hugo_symbol == "STK11") %>%
  mutate(num_stk11_samples = n_distinct(tumor_sample_barcode)) %>%
  group_by(group) %>%
  mutate(
    group_num_stk11_samples = n_distinct(tumor_sample_barcode),
    group_num_stk11_mutations = n()
  ) %>%
  group_by(
    group, num_luad_samples, num_stk11_samples,
    amino_position, mutation_type,
    group_num_stk11_samples, group_num_stk11_mutations
  ) %>%
  summarise(
    group_position_type_num_stk11_samples = n_distinct(tumor_sample_barcode),
    num_mutation_type = n(),
    percent_of_mutations = num_mutation_type / unique(group_num_stk11_mutations),
    sifts = paste(unique(sift_pred), collapse = ", "),
    polyphens = paste(unique(polyphen2_hdiv_pred), collapse = ", "),
    fathmms = paste(unique(fathmm_pred), collapse = ", "),
    clinsigs = paste(unique(clinsig), collapse = ", ")
  ) %>%
  mutate(
    sifts = parse_sift(sifts),
    polyphens = parse_polyphen(polyphens),
    fathmms = parse_fathmm(fathmms),
    clinsigs = parse_clinsig(clinsigs)
  ) %>%
  ungroup()

stk11_mutations %<>%
  mutate(mutation_type = ifelse(mutation_type == "frameshift_deletion",
    "frame_shift_del", mutation_type
  ))

stk11_mutations_dmg <- stk11_mutations %>%
  group_by(group, amino_position) %>%
  summarise(
    num_mutations = sum(num_mutation_type),
    damaging = any(sifts, polyphens, fathmms, clinsigs)
  ) %>%
  ungroup() %>%
  filter(damaging) %>%
  mutate(
    amino_position = as.numeric(amino_position),
    num_mutations = ifelse(group == "G12C",
      num_mutations,
      -num_mutations
    )
  )


stk11_lollipop <- stk11_mutations %>%
  mutate(
    amino_position = as.numeric(amino_position),
    num_mutation_type = ifelse(group == "G12C",
      num_mutation_type,
      -num_mutation_type
    )
  ) %>%
  ggplot(aes(x = amino_position, y = num_mutation_type)) +
  geom_col(aes(fill = mutation_type), position = "stack") +
  geom_point(
    data = stk11_mutations_dmg,
    aes(x = amino_position, y = num_mutations),
    color = "tomato", size = 1
  ) +
  theme_bw(base_size = 7, base_family = "Arial")
ggsave_wrapper(
  stk11_lollipop,
  plot_path(GRAPHS_DIR, "stk11_lollipop.svg"),
  "wide"
)



#### ---- Statistics ---- ####

# Fisher test for difference in types of mutations.
stk11_mutations %>%
  filter(!is.na(amino_position)) %>%
  select(group, amino_position, mutation_type) %>%
  count(group, mutation_type, name = "mut_type_total") %>%
  pivot_wider(c(group, mutation_type, mut_type_total),
    names_from = group,
    values_from = mut_type_total
  ) %>%
  as.data.frame() %>%
  column_to_rownames("mutation_type") %>%
  fisher.test()

# >     Fisher's Exact Test for Count Data
# >
# > data:  .
# > p-value = 0.726
# > alternative hypothesis: two.sided


#### ---- STK11 Lollipop ---- ####

add_damage_pts <- function(p, group,
                           color = "firebrick", size = 0.6, alpha = 1.0) {
  dmg_df <- stk11_mutations_dmg %>%
    filter(group == !!group) %>%
    filter(!is.na(amino_position))
  p <- p + geom_point(
    data = dmg_df,
    aes(x = amino_position, y = num_mutations),
    color = color, size = size, alpha = alpha
  )
  return(p)
}


stk11_pal <- c(
  "N-term" = "goldenrod1",
  "NLS" = "mediumorchid1",
  "Kinase" = "darkturquoise",
  "C-term" = "darksalmon"
)

stk11_scheme <- tibble(
  x = c(1, 38, 38, 43, 43, 49, 49, 309, 309, 433),
  domain = c(
    "N-term", "N-term", "NLS", "NLS", "N-term", "N-term",
    "Kinase", "Kinase", "C-term", "C-term"
  )
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
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
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


stk11_breaks <- c(1, seq(25, 401, 25), 433)

stk11_g12c_lollipop <- stk11_mutations %>%
  filter(group == "G12C") %>%
  mutate(
    amino_position = as.numeric(amino_position),
    mutation_type = format_var_names(mutation_type)
  ) %>%
  ggplot(aes(x = amino_position, y = num_mutation_type)) +
  geom_col(aes(fill = mutation_type), position = "stack") +
  geom_hline(yintercept = 0, size = 0.8, color = "black") +
  scale_fill_manual(values = mod_variant_pal) +
  scale_x_continuous(
    limits = c(min(stk11_breaks), max(stk11_breaks)),
    expand = expansion(mult = c(0.01, 0.01)),
    breaks = stk11_breaks
  ) +
  scale_y_continuous(expand = expansion(add = c(0, 0.15))) +
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
stk11_g12c_lollipop <- add_damage_pts(stk11_g12c_lollipop, "G12C")
ggsave_wrapper(
  stk11_g12c_lollipop,
  plot_path(GRAPHS_DIR, "stk11_g12c_lollipop.svg"),
  "small"
)

stk11_rest_lollipop <- stk11_mutations %>%
  filter(group != "G12C") %>%
  mutate(
    amino_position = as.numeric(amino_position),
    mutation_type = format_var_names(mutation_type),
    num_mutation_type = -1 * num_mutation_type
  ) %>%
  ggplot(aes(x = amino_position, y = num_mutation_type)) +
  geom_col(aes(fill = mutation_type), position = "stack") +
  geom_hline(yintercept = 0, size = 0.8, color = "black") +
  scale_fill_manual(values = mod_variant_pal) +
  scale_x_continuous(
    limits = c(min(stk11_breaks), max(stk11_breaks)),
    expand = expansion(mult = c(0.01, 0.01)),
    breaks = stk11_breaks
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.15, 0)),
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
    x = "STK11 amino acid position"
  )
stk11_rest_lollipop <- add_damage_pts(stk11_rest_lollipop, "rest")
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

saveFigRds(stk11_lollipop_patch, "stk11_lollipop_patch.rds")
