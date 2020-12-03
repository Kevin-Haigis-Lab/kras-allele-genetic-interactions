# Plot the distribution of KRAS mutations.

GRAPHS_DIR <- "90_05_kras-allele-distribution"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)

# Number of WES/WGS tumor samples per cancer.
cancer_full_muts_df %>%
  filter(target %in% c("genome", "exome")) %>%
  distinct(cancer, tumor_sample_barcode) %>%
  count(cancer) %>%
  knitr::kable(format = "simple")
# > cancer       n
# > -------  -----
# > COAD      1536
# > LUAD       891
# > MM        1201
# > PAAD      1395
# > SKCM      1042

# A data frame of the KRAS allele for each `tumor_sample_barcode`.
alleles_df <- cancer_full_muts_df %>%
  group_by(cancer, tumor_sample_barcode) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    cancer, tumor_sample_barcode, dataset, target,
    is_hypermutant, ras, ras_allele
  ) %>%
  unique()

# The alleles to use for plotting.
alleles_to_keep <- names(short_allele_pal)
alleles_to_keep <- alleles_to_keep[alleles_to_keep != "Other"]

# A data frame of the distribution of alleles across cancer.
get_allele_distribtion_dataframe <- function(with_other = TRUE) {
  df <- alleles_df %>%
    filter(!is_hypermutant) %>%
    mutate(ras_allele = str_remove(ras_allele, "KRAS_"))

  if (with_other) {
    df %<>%
      mutate(
        ras_allele = fct_other(ras_allele, keep = alleles_to_keep)
      ) %>%
      group_by(cancer) %>%
      mutate(
        ras_allele = as.character(fct_lump(ras_allele, prop = 0.01)),
      )
  }
  df %<>%
    group_by(cancer) %>%
    mutate(
      num_cancer_samples = n_distinct(tumor_sample_barcode)
    ) %>%
    group_by(cancer, ras, ras_allele, num_cancer_samples) %>%
    summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
    ungroup()

  df %<>%
    filter(cancer != "SKCM") %>%
    mutate(allele_freq = num_allele_samples / num_cancer_samples)

  return(df)
}

allele_dist <- get_allele_distribtion_dataframe(with_other = TRUE)
allele_dist_all <- get_allele_distribtion_dataframe(FALSE)


# Return the alleles in order of the frequency, except with "Other" always
# at the end.
get_factor_levels <- function(allele, freq) {
  lvls <- allele[order(-freq)]
  if (any(allele == "Other")) {
    lvls <- c(lvls[lvls != "Other"], "Other")
  }
  return(lvls)
}


# Make a bar plot for the distribution of the alleles.
# Provide a value to `max_freq` to set the y-axis limit.
make_allele_dist_barplot <- function(cancer, data,
                                     max_freq = NA,
                                     lvls = NULL) {
  data <- data %>% filter(ras_allele != "WT")

  factor_levels <- unique(get_factor_levels(data$ras_allele, data$allele_freq))
  if (!is.null(lvls)) factor_levels <- unique(lvls)

  p <- data %>%
    mutate(ras_allele = factor(ras_allele, levels = !!factor_levels)) %>%
    ggplot(aes(x = ras_allele, y = allele_freq)) +
    geom_col(aes(fill = ras_allele)) +
    scale_fill_manual(
      values = short_allele_pal,
      guide = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, max_freq),
      expand = expansion(mult = c(0, 0.02)),
      breaks = round_breaks()
    ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank()
    ) +
    coord_flip() +
    labs(
      title = cancer,
      y = "frequency"
    )

  return(p)
}


# Make a stacked bar plot of the allele frequency.
make_allele_stackedplot <- function(cancer, data, ...) {
  factor_levels <- get_factor_levels(data$ras_allele, data$allele_freq)
  factor_levels <- c("WT", factor_levels[factor_levels != "WT"])
  p <- data %>%
    mutate(ras_allele = factor(ras_allele, levels = !!factor_levels)) %>%
    ggplot(aes(x = "KRAS", y = allele_freq)) +
    geom_col(aes(fill = ras_allele), position = "stack") +
    scale_fill_manual(
      values = short_allele_pal,
      guide = FALSE
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 1),
      breaks = seq(0, 1.0, 0.25)
    ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank()
    )

  return(p)
}


# A wrapper to save the bar plot for a cancer.
save_allele_dist_barplot <- function(barplot, cancer,
                                     size = NA, width = NA, height = NA,
                                     ...) {
  ggsave_wrapper(
    barplot,
    plot_path(
      GRAPHS_DIR,
      glue("allele_dist_barplot_{cancer}.svg")
    ),
    size, width, height
  )
}


max_freq <- allele_dist %>%
  filter(ras_allele != "WT") %>%
  pull(allele_freq) %>%
  max()

plots <- allele_dist %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(
    barplot = purrr::map2(cancer, data, make_allele_dist_barplot,
      max_freq = NA
    ),
    stackedplot = purrr::map2(cancer, data, make_allele_stackedplot)
  ) %>%
  pwalk(save_allele_dist_barplot, width = 4, height = 2.5)

ggsave_wrapper(
  wrap_plots(plots$barplot),
  plot_path(
    GRAPHS_DIR,
    glue("allele_dist_barplot_all.svg")
  ),
  width = 8, height = 4.5
)


barplots_facet <- make_allele_dist_barplot(
  cancer = "all",
  data = allele_dist,
  lvls = names(short_allele_pal)
) +
  facet_wrap(~cancer, scales = "free_x", nrow = 2) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02)),
    breaks = round_breaks(),
    labels = function(x) ifelse(x == 0.0, "", x)
  ) +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9)
  ) +
  labs(
    y = "fraction of samples",
    x = "*KRAS* alleles"
  )
ggsave_wrapper(
  barplots_facet,
  plot_path(GRAPHS_DIR, glue("allele_dist_barplot_all_facet.jpeg")),
  width = 6.5, height = 3.5
)


plot_distribution_and_stacked <- function(cancer, data, with_extra_space = FALSE) {
  bp <- make_allele_dist_barplot(cancer, data)
  sp <- make_allele_stackedplot(cancer, data)

  p <- bp + sp + plot_layout(widths = c(20, 1))

  if (with_extra_space) {
    p <- p + plot_spacer() + plot_layout(widths = c(20, 1, 1))
  }

  return(p)
}

plots <- allele_dist %>%
  group_by(cancer) %>%
  nest() %>%
  ungroup() %>%
  arrange(cancer) %>%
  mutate(with_extra_space = rep(c(TRUE, FALSE), 2)) %>%
  pmap(plot_distribution_and_stacked)
ggsave_wrapper(
  patchwork::wrap_plots(plots),
  plot_path(
    GRAPHS_DIR,
    glue("allele_dist_barplot_stackplot.svg")
  ),
  width = 8, height = 4.5
)

saveFigRds(plots, "allele_dist_barplot_stackplot")


# Plot all of the alleles (with at least 3 appearences)
# Need to add the new colors to the palette - I just changed the
#   `short_allele_pal`, though, if this analysis gets any bigger, then I will
#   parameterize the palette for the plotting functions (or pass in `...`).
{
  ORIGINAL_PAL <- short_allele_pal
  added_alleles <- setdiff(unique(allele_dist_all$ras_allele), names(short_allele_pal))
  short_allele_pal <- c(
    short_allele_pal,
    rep(short_allele_pal["Other"], length(added_alleles))
  )
  names(short_allele_pal) <- c(names(ORIGINAL_PAL), added_alleles)

  plots <- allele_dist_all %>%
    filter(num_allele_samples > 2) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    arrange(cancer) %>%
    mutate(with_extra_space = rep(c(TRUE, FALSE), 2)) %>%
    pmap(plot_distribution_and_stacked)
  ggsave_wrapper(
    patchwork::wrap_plots(plots),
    plot_path(
      GRAPHS_DIR,
      glue("allele_dist_barplot_stackplot_all.svg")
    ),
    width = 8, height = 4.5
  )

  saveFigRds(plots, "allele_dist_barplot_stackplot_all")

  short_allele_pal <- ORIGINAL_PAL
}


#### ---- Dot-plot of allele frequency ---- ####


cancer_to_longname <- tibble(
  cancer = c("COAD", "LUAD", "MM", "PAAD", "SKCM"),
  long_cancer = c(
    "COAD\n(colon", "LUAD\n(lung", "MM\n(WBC",
    "PAAD\n(pancreas", "SKCM\n(skin"
  )
) %>%
  left_join(distinct(select(allele_dist, cancer, num_cancer_samples)),
    by = "cancer"
  ) %>%
  mutate(
    num = scales::comma(num_cancer_samples),
    long_cancer = paste0(long_cancer, ", n=", num, ")"),
    long_cancer = fct_rev(fct_inorder(long_cancer))
  ) %>%
  select(cancer, long_cancer)

long_cancer_palette <- cancer_palette
names(long_cancer_palette) <- cancer_to_longname$long_cancer


allele_dist_dotplot <- allele_dist_all %>%
  filter(allele_freq > 0.001 & num_allele_samples >= 5) %>%
  filter(!(ras_allele %in% c("WT", "Other", "K117N")) &
    ras_allele %in% names(short_allele_pal)) %>%
  group_by(cancer) %>%
  mutate(
    allele_freq = num_allele_samples / sum(num_allele_samples),
    codon = str_extract(ras_allele, "[:digit:]+"),
    codon = as.numeric(codon)
  ) %>%
  ungroup() %>%
  left_join(cancer_to_longname, by = "cancer") %>%
  ggplot(aes(x = ras_allele, y = long_cancer)) +
  facet_grid(. ~ codon, scales = "free_x", space = "free") +
  geom_point(aes(size = allele_freq, color = long_cancer)) +
  scale_size_area() +
  scale_color_manual(values = long_cancer_palette, guide = NULL) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_markdown(),
    axis.title.y = element_blank(),
    legend.position = "left",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "*KRAS* alleles")

ggsave_wrapper(
  allele_dist_dotplot,
  plot_path(GRAPHS_DIR, "allele_dist_dotplot.svg"),
  "wide"
)

saveFigRds(allele_dist_dotplot, "allele_dist_dotplot")


#### ---- Barplot of KRAS mut freq per cancer ---- ####

cancer_freq_kras_mut <- cancer_full_muts_df %>%
  filter(!is_hypermutant & cancer != "SKCM") %>%
  group_by(cancer, tumor_sample_barcode, ras_allele) %>%
  slice(1) %>%
  ungroup() %>%
  select(cancer, tumor_sample_barcode, ras_allele) %>%
  group_by(cancer) %>%
  summarise(freq_kras_mut = sum(ras_allele != "WT") / n()) %>%
  ungroup() %>%
  mutate(cancer = factor(cancer, levels = rev(names(cancer_palette))))
knitr::kable(cancer_freq_kras_mut, digits = 3, format = "simple")
# > cancer    freq_kras_mut
# > -------  --------------
# > COAD              0.414
# > LUAD              0.353
# > MM                0.219
# > PAAD              0.863

cancer_freq_kras_mut_column <- cancer_freq_kras_mut %>%
  left_join(cancer_to_longname, by = "cancer") %>%
  ggplot(aes(x = freq_kras_mut, y = long_cancer)) +
  geom_col(aes(fill = cancer), width = 0.6) +
  scale_x_continuous(expand = expansion(add = c(0, 0.02))) +
  scale_fill_manual(values = cancer_palette) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_markdown(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "freq. of *KRAS* mut.")
ggsave_wrapper(
  cancer_freq_kras_mut_column,
  plot_path(GRAPHS_DIR, "cancer_freq_kras_mut_column.svg"),
  "small"
)
saveFigRds(cancer_freq_kras_mut_column, "cancer_freq_kras_mut_column")


#### ---- Lollipop plot of KRAS mutations ---- ####

codons_to_label <- c(12, 13, 61, 146)

# 'maftools' lollipop plot. (only if X11 is available)
caps <- as.list(capabilities())
if (caps$X11) {
  kras_maf <- cancer_full_coding_muts_maf %>%
    filter(hugo_symbol == "KRAS") %>%
    maftools::read.maf(verbose = FALSE)

  svg(plot_path(GRAPHS_DIR, "lollipop-kras.svg"), width = 6, height = 4)
  maftools::lollipopPlot(kras_maf,
    "KRAS",
    AACol = "amino_position",
    labelPos = codons_to_label,
    titleSize = c(0.1, 0.1),
    pointSize = 1.2,
    axisTextSize = c(0.75, 0.75),
    showDomainLabel = FALSE
  )
  dev.off()
}

# My lollipop plot.
kras_lollipop_plot <- cancer_full_coding_muts_maf %>%
  filter(!is_hypermutant & hugo_symbol == "KRAS") %>%
  mutate(amino_position = as.numeric(amino_position)) %>%
  filter(!is.na(amino_position)) %>%
  group_by(cancer, amino_position) %>%
  summarise(num_apos = n_distinct(tumor_sample_barcode)) %>%
  group_by(amino_position) %>%
  mutate(
    total_num_apos = sum(num_apos),
    log_total_num_apos = log10(total_num_apos)
  ) %>%
  ungroup() %>%
  mutate(
    cancer_frac_apos = num_apos / total_num_apos,
    cancer_frac_log_apos = cancer_frac_apos * log_total_num_apos,
    point_label = ifelse(
      amino_position %in% !!codons_to_label & cancer == "COAD",
      as.character(amino_position),
      NA
    )
  ) %>%
  ggplot(aes(x = amino_position)) +
  geom_col(
    aes(
      y = cancer_frac_log_apos,
      fill = cancer
    )
  ) +
  geom_point(
    aes(y = log_total_num_apos),
    color = "grey25",
    size = 0.8
  ) +
  geom_text(
    aes(
      label = point_label,
      y = log_total_num_apos
    ),
    family = "Arial",
    size = 2,
    hjust = 0,
    nudge_x = 5,
    nudge_y = 0
  ) +
  scale_fill_manual(
    values = cancer_palette,
    guide = guide_legend(
      title = NULL,
      label.hjust = 0,
      keywidth = unit(2, "mm"),
      keyheight = unit(2, "mm"),
      ncol = 1
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 4),
    breaks = c(0, 1, 2, 3, 4, 5),
    labels = function(x) {
      10^x
    }
  ) +
  scale_x_continuous(
    limits = c(1, 189),
    expand = c(0, 0)
  ) +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = c(0.8, 0.8),
    legend.spacing.x = unit(1, "mm"),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    axis.title.y = element_markdown()
  ) +
  labs(
    x = "KRas amino acid sequence",
    y = "number of samples (*log*<sub>10</sub>-transformed)"
  )

ggsave_wrapper(
  kras_lollipop_plot,
  plot_path(GRAPHS_DIR, "lollipop-kras_2.svg"),
  "small"
)

# Save for use in Figure 1.
saveFigRds(kras_lollipop_plot, "lollipop-kras_2")


#### ---- Lollipop plot of KRAS mutations (only hot-spot) ---- ####

manual_spacing_levels <- c(1, 12, 13, 15:16, 61, 65:67, 146, 188:189)

# My lollipop plot.
kras_lollipop_plot <- cancer_full_coding_muts_maf %>%
  filter(!is_hypermutant & hugo_symbol == "KRAS") %>%
  mutate(amino_position = as.numeric(amino_position)) %>%
  filter(!is.na(amino_position)) %>%
  filter(amino_position %in% codons_to_label) %>%
  group_by(cancer, amino_position) %>%
  summarise(num_apos = n_distinct(tumor_sample_barcode)) %>%
  group_by(amino_position) %>%
  mutate(
    total_num_apos = sum(num_apos),
    log_total_num_apos = log10(total_num_apos)
  ) %>%
  ungroup() %>%
  mutate(
    cancer_frac_apos = num_apos / total_num_apos,
    cancer_frac_log_apos = cancer_frac_apos * log_total_num_apos,
    point_label = ifelse(
      amino_position %in% !!codons_to_label & cancer == "COAD",
      as.character(amino_position),
      NA
    ),
    amino_position = factor(amino_position, levels = manual_spacing_levels)
  ) %>%
  ggplot(aes(x = amino_position)) +
  geom_col(
    aes(
      y = cancer_frac_log_apos,
      fill = cancer
    )
  ) +
  scale_fill_manual(
    values = cancer_palette,
    guide = guide_legend(
      title = NULL,
      label.hjust = 0,
      keywidth = unit(3, "mm"),
      keyheight = unit(3, "mm"),
      ncol = 1
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 4),
    breaks = c(0, 1, 2, 3, 4, 5),
    labels = function(x) {
      10^x
    }
  ) +
  scale_x_discrete(
    expand = c(0, 0), drop = FALSE,
    labels = function(x) {
      ifelse(x %in% c(codons_to_label, 1, 189), x, "")
    }
  ) +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = c(0.8, 0.8),
    legend.spacing.x = unit(1, "mm"),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    axis.title.y = element_markdown(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "KRas amino acid sequence",
    y = "number of samples (*log*<sub>10</sub>-transformed)"
  )

ggsave_wrapper(
  kras_lollipop_plot,
  plot_path(GRAPHS_DIR, "lollipop-kras_hotspot-only.svg"),
  "small"
)

# Save for use in Figure 1.
saveFigRds(kras_lollipop_plot, "lollipop-kras_hotspot-only")



#### ---- Table of the distribution of alleles ---- ####

alleles_df %<>%
  filter(cancer != "SKCM" & !is_hypermutant) %>%
  mutate(
    kras_allele = str_remove(ras_allele, "KRAS_"),
    codon = str_extract(ras_allele, "[:digit:]+|WT")
  )

# Calculate the frequency of the alleles of KRAS by cancer.
calc_frequency_of_alleles_by_cancer <- function(df) {
  df %>%
    group_by(cancer) %>%
    mutate(
      num_cancer_samples = n_distinct(tumor_sample_barcode)
    ) %>%
    group_by(cancer, ras, kras_allele, num_cancer_samples) %>%
    summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    mutate(allele_frequency = num_allele_samples / num_cancer_samples) %>%
    arrange(cancer, -allele_frequency) %>%
    select(
      cancer, kras_allele,
      num_allele_samples, num_cancer_samples, allele_frequency
    )
}

# Frequency of each allele across cancers.
alleles_df %>%
  calc_frequency_of_alleles_by_cancer() %T>%
  write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution.tsv")) %>%
  save_supp_data(3, "KRAS allele frequencies")

# Frequency of each allele across cancers without WT.
alleles_df %>%
  filter(kras_allele != "WT") %>%
  calc_frequency_of_alleles_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution-noWT.tsv"))


# Frequency of each allele for all cancers without WT.
alleles_df %>%
  mutate(cancer = "ALL") %>%
  calc_frequency_of_alleles_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution-all.tsv"))

# Frequency of each allele for all cancers without WT.
alleles_df %>%
  mutate(cancer = "ALL") %>%
  filter(kras_allele != "WT") %>%
  calc_frequency_of_alleles_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution-all-noWT.tsv"))


# Calculate the frequency of mutations of KRAS codons by cancer.
calc_frequency_of_codons_by_cancer <- function(df) {
  df %>%
    group_by(cancer) %>%
    mutate(
      num_cancer_samples = n_distinct(tumor_sample_barcode)
    ) %>%
    group_by(cancer, codon, num_cancer_samples) %>%
    summarise(num_codon_samples = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    mutate(codon_frequency = num_codon_samples / num_cancer_samples) %>%
    arrange(cancer, -codon_frequency) %>%
    select(
      cancer, codon,
      num_codon_samples, num_cancer_samples, codon_frequency
    )
}

# The frequency of each codon across cancers.
alleles_df %>%
  calc_frequency_of_codons_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution.tsv"))

# The frequency of each codon across cancers without WT.
alleles_df %>%
  filter(kras_allele != "WT") %>%
  calc_frequency_of_codons_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution-noWT.tsv"))


# The frequency of codon with all cancers combined.
alleles_df %>%
  mutate(cancer = "ALL") %>%
  calc_frequency_of_codons_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution-all.tsv"))

# The frequency of codon with all cancers combined without WT.
alleles_df %>%
  filter(kras_allele != "WT") %>%
  mutate(cancer = "ALL") %>%
  calc_frequency_of_codons_by_cancer() %>%
  write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution-all-noWT.tsv"))

# The frequency of codon with all cancers combined without WT and K117.
alleles_df %>%
  filter(!str_detect(kras_allele, "117|WT")) %>%
  mutate(cancer = "ALL") %>%
  calc_frequency_of_codons_by_cancer() %>%
  write_tsv(
    table_path(GRAPHS_DIR, "kras-codon-distribution-all-noWTK117.tsv")
  )



#### ---- Accounting for cancer prevelance ---- ####

incidence_source_pal <- c(
  "ACS" = "#B60A2D",
  "SEER" = "#0A57A6"
)

cancer_incidence_plot <- cancer_incidence_df %>%
  filter(cancer != "SKCM") %>%
  ggplot(aes(x = cancer, y = incidence)) +
  geom_col(aes(fill = source), position = "dodge") +
  scale_fill_manual(values = incidence_source_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "grey50")
  ) +
  labs(
    y = "average yearly incidence (2012-2016)"
  )
ggsave_wrapper(
  cancer_incidence_plot,
  plot_path(GRAPHS_DIR, "cancer_incidence_plot.svg"),
  "small"
)


cancer_incidence_stacked_plot <- cancer_incidence_df %>%
  filter(cancer != "SKCM") %>%
  filter(cancer != "all") %>%
  ggplot(aes(x = source, y = incidence)) +
  geom_col(aes(fill = cancer), position = "stack") +
  scale_fill_manual(values = cancer_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
  ) +
  labs(
    y = "average yearly incidence (2012-2016)"
  )
ggsave_wrapper(
  cancer_incidence_stacked_plot,
  plot_path(GRAPHS_DIR, "cancer_incidence_stacked_plot.svg"),
  "small"
)


cancer_incidence_diff_plot <- cancer_incidence_df %>%
  filter(cancer != "SKCM") %>%
  pivot_wider(
    id_cols = cancer, names_from = source, values_from = incidence
  ) %>%
  mutate(
    source_diff = ACS - SEER,
    fill_val = ifelse(source_diff < 0, -1, 1) * log10(abs(source_diff))
  ) %>%
  ggplot(aes(x = cancer, y = source_diff)) +
  geom_col(aes(fill = fill_val)) +
  scale_fill_gradient2(
    low = "dodgerblue", mid = "grey95", high = "tomato",
    na.value = "white"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  labs(
    y = "difference in incidence, ACS - SEER"
  )
ggsave_wrapper(
  cancer_incidence_diff_plot,
  plot_path(GRAPHS_DIR, "cancer_incidence_diff_plot.svg"),
  "small"
)


percent_of_cancer <- tibble::enframe(FRACTION_OF_TISSUE_THAT_ARE_CANCER,
  name = "cancer",
  value = "frac"
) %>%
  unnest(frac)

# The incidence weights for each cancer.
incidence_data <- acs_incidence_df %>%
  filter(!is.na(cancer) & !(cancer %in% c("SKCM", "all"))) %>%
  select(cancer, incidence) %>%
  left_join(percent_of_cancer, by = "cancer") %>%
  mutate(
    incidence = incidence * frac,
    incidence_frac_of_cancers = incidence / sum(incidence)
  ) %>%
  select(-frac)

# Fraction of KRAS mutations at the hotspot codons.
# `avg_codon_freq` is the average of the frequencies across the 4 cancers
# `adj_codon_freq` is the average, weighted by prevalence of the cancer
alleles_df %>%
  filter(kras_allele != "WT") %>%
  calc_frequency_of_codons_by_cancer() %>%
  left_join(incidence_data, by = "cancer") %>%
  select(
    cancer, codon, codon_frequency,
    incidence, incidence_frac_of_cancers
  ) %>%
  mutate(
    cancer_codon_freq = codon_frequency * incidence_frac_of_cancers,
    cancer_codon_freq = softmax(cancer_codon_freq)
  ) %>%
  select(cancer, codon, codon_frequency, cancer_codon_freq) %>%
  group_by(codon) %>%
  summarise(
    avg_codon_freq = sum(codon_frequency),
    adj_codon_freq = sum(cancer_codon_freq)
  ) %>%
  ungroup() %>%
  filter(codon != "117") %>%
  mutate(
    avg_codon_freq = softmax(avg_codon_freq),
    adj_codon_freq = softmax(adj_codon_freq)
  ) %>%
  arrange(-adj_codon_freq) %T>%
  write_tsv(table_path(GRAPHS_DIR, "fraction-kras-percodon-adjusted.tsv")) %>%
  knitr::kable(digits = 3)
# > |codon | avg_codon_freq| adj_codon_freq|
# > |:-----|--------------:|--------------:|
# > |12    |          0.722|          0.768|
# > |13    |          0.098|          0.114|
# > |61    |          0.148|          0.081|
# > |146   |          0.032|          0.037|


#### ---- Association of alleles with hypermutation status ---- ####
# This analysis is only relevant in COAD.

kras_hypermuts <- cancer_full_coding_muts_df %>%
  filter(cancer == "COAD") %>%
  distinct(tumor_sample_barcode, is_hypermutant, ras_allele)
kras_hypermuts


make_lgl_table <- function(x, y) {
  f <- function(a) {
    factor(a, levels = c("FALSE", "TRUE"))
  }
  table(f(x), f(y))
}


test_for_hypermutant_association <- function(allele, data) {
  data %>%
    mutate(is_allele = ras_allele == !!allele) %$%
    make_lgl_table(is_allele, is_hypermutant) %>%
    fisher.test(alternative = "g") %>%
    broom::tidy()
}


kras_alleles_hypermut_association <- kras_hypermuts %>%
  count(ras_allele) %>%
  filter(n > 20) %>%
  select(-n) %>%
  mutate(hypermut_ass_res = map(ras_allele,
    test_for_hypermutant_association,
    data = kras_hypermuts
  ))
kras_alleles_hypermut_association %>%
  unnest(hypermut_ass_res) %>%
  janitor::clean_names() %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  arrange(ras_allele) %>%
  select(-method, -alternative) %>%
  mutate(
    log_or = log(estimate),
    sig = ifelse(estimate > 1 & adj_p_value < 0.05, "X", " ")
  ) %>%
  knitr::kable(digits = 3, format = "simple")
# > ras_allele    estimate   p_value   conf_low   conf_high   adj_p_value   log_or  sig
# > -----------  ---------  --------  ---------  ----------  ------------  -------  ----
# > KRAS_A146T       1.599     0.024      1.082         Inf         0.071    0.469
# > KRAS_A146V       1.307     0.372      0.471         Inf         0.744    0.267
# > KRAS_G12A        0.944     0.622      0.535         Inf         0.933   -0.058
# > KRAS_G12C        0.238     1.000      0.092         Inf         1.000   -1.436
# > KRAS_G12D        1.009     0.490      0.816         Inf         0.840    0.009
# > KRAS_G12S        0.531     0.970      0.241         Inf         1.000   -0.633
# > KRAS_G12V        0.319     1.000      0.215         Inf         1.000   -1.141
# > KRAS_G13D        1.269     0.059      0.987         Inf         0.143    0.238
# > KRAS_K117N       3.219     0.005      1.518         Inf         0.030    1.169  X
# > KRAS_Q61H        0.178     0.995      0.009         Inf         1.000   -1.724
# > KRAS_Q61K        6.347     0.000      2.922         Inf         0.000    1.848  X
# > WT               1.219     0.010      1.059         Inf         0.038    0.198  X

kras_hypermuts %>%
  count(ras_allele, sort = TRUE) %>%
  mutate(freq = n / sum(n)) %>%
  knitr::kable(digits = 3, format = "simple")
# > ras_allele        n    freq
# > ------------  -----  ------
# > WT             2872   0.592
# > KRAS_G12D       589   0.121
# > KRAS_G12V       410   0.084
# > KRAS_G13D       359   0.074
# > KRAS_A146T      130   0.027
# > KRAS_G12C       123   0.025
# > KRAS_G12A        92   0.019
# > KRAS_G12S        82   0.017
# > KRAS_Q61H        33   0.007
# > KRAS_K117N       28   0.006
# > KRAS_A146V       27   0.006
# > KRAS_Q61K        23   0.005
# > KRAS_G12R        20   0.004
# > KRAS_Q61L        19   0.004
# > KRAS_G13C        17   0.004
# > KRAS_Q61R        17   0.004
# > KRAS_G13R         2   0.000
# > KRAS_K117R        2   0.000
# > KRAS_Q61E         2   0.000
# > KRAS_Q61P         2   0.000
# > KRAS_A146P        1   0.000
# > KRAS_G12F         1   0.000
# > KRAS_G13dup       1   0.000
# > KRAS_GV13CG       1   0.000
