library(mustashe)
library(jhcutils)
library(glue)
library(ggridges)
library(patchwork)
library(ggtext)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
  if (str_detect(f, "global|enrich")) {
    next
  }
  source(f)
}

options(dplyr.summarise.inform = FALSE)


kras_allele_causation_mutsig_df <- readRDS("~/Downloads/kras_allele_causation_mutsig_df.rds")
mutsig_noartifact_df <- readRDS("~/Downloads/mutsig_noartifact_df.rds")
allele_signature_associations <- readRDS("~/Downloads/allele_signature_associations.rds")

theme_set(
  theme_bw(base_size = 11, base_family = "Arial") %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.ticks = element_blank()
    )
)



alleles_tested_per_cancer <- allele_signature_associations %>%
  group_by(cancer) %>%
  summarise(
    alleles_tested = list(unique(c(unlist(group1), unlist(group2))))
  ) %>%
  ungroup()



ggridge_plot_signatures_pgridge <- function(ms_df, signature) {
  ms_df %>%
    ggplot(aes(x = contribution, y = allele)) +
    geom_density_ridges(
      aes(color = allele, fill = allele),
      alpha = 0.4
    ) +
    scale_color_manual(
      values = short_allele_pal,
      guide = FALSE
    ) +
    scale_fill_manual(
      values = short_allele_pal,
      guide = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0.02))
    ) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      axis.text.y = element_text(vjust = 0)
    ) +
    labs(
      x = glue("mutational signature {signature} level"),
      y = NULL
    )
}


orderless_between <- function(x, a, b) {
  hold_a <- a
  a <- if_else(a <= b, a, b)
  b <- if_else(hold_a < b, b, hold_a)
  map2_lgl(a, b, function(a, b) {
    between(x, a, b)
  })
}


value_overlaps_existing_placement <- function(stats_df, idx, current_x) {
  idx_g1 <- stats_df$g1[stats_df$rank == idx]
  idx_g2 <- stats_df$g2[stats_df$rank == idx]
  overlap_df <- stats_df %>%
    filter(rank < !!idx & x == !!current_x)

  if (nrow(overlap_df) == 0) {
    return(FALSE)
  }

  overlap_df <- overlap_df %>%
    mutate(is_overlap = orderless_between(idx_g1, g1, g2) | orderless_between(idx_g2, g1, g2)) %>%
    filter(is_overlap)

  return(nrow(overlap_df) != 0)
}


optimize_stats_bar_placement <- function(stats_df) {
  stats_df$x <- 0
  for (i in seq(1, nrow(stats_df))) {
    if (i == 1) {
      next
    }
    x <- 0
    is_overlap <- TRUE
    while (is_overlap & x < nrow(stats_df)) {
      if (!value_overlaps_existing_placement(stats_df, i, x)) {
        break
      }
      x <- x + 1
    }
    stats_df$x[[i]] <- x
  }
  stats_df
}


ggridge_plot_signatures_pstats <- function(stats_df,
                                           alleles,
                                           added_height = 0,
                                           extra_x_left = 0.1) {
  stats_df$id <- seq(1, nrow(stats_df))

  plot_df <- stats_df %>%
    mutate(
      group1 = factor(group1, levels = alleles),
      group2 = factor(group2, levels = alleles),
      g1 = as.numeric(group1),
      g2 = as.numeric(group2),
      lowest_g = map2_dbl(g1, g2, ~ min(c(.x, .y))),
      dist_g = abs(g1 - g2)
    ) %>%
    arrange(lowest_g, dist_g) %>%
    mutate(rank = row_number()) %>%
    optimize_stats_bar_placement() %>%
    mutate(
      p_value_star = assign_stars(p_value),
      star_x = x,
      star_y = map2_dbl(g1, g2, ~ mean(c(.x, .y)))
    ) %>%
    mutate(
      x = -1 * x,
      star_x = -1 * star_x
    )

  x_lims <- c(
    min(plot_df$star_x) - extra_x_left,
    max(plot_df$x)
  )

  y_lims <- c(
    1,
    max(c(plot_df$g1, plot_df$g2)) + added_height
  )

  plot_df %>%
    select(x, g1, g2, id) %>%
    pivot_longer(-c(id, x), names_to = NULL, values_to = "group_value") %>%
    ggplot(aes(x = x, y = group_value)) +
    geom_line(aes(group = id)) +
    geom_text(
      aes(x = star_x, y = star_y, label = p_value_star),
      data = plot_df,
      angle = 90,
      hjust = 0.5,
      vjust = 0,
    ) +
    scale_x_continuous(
      limits = x_lims,
      expand = expansion(mult = c(0.15, 0.02))
    ) +
    scale_y_continuous(
      limits = y_lims,
      expand = expansion(mult = c(0, 0.02)),
      breaks = seq(0, y_lims[[2]])
    ) +
    theme_void(base_size = 7, base_family = "Arial")
}



ggridge_plot_signatures <- function(signature,
                                    cancer,
                                    stats_df,
                                    ms_df,
                                    alleles_tested = NULL,
                                    stats_plot_add_height = 1.5,
                                    patch_widths = c(1, 10)) {
  mod_ms_df <- ms_df %>%
    filter(cancer == !!cancer & signature == !!signature) %>%
    mutate(
      allele = ifelse(ras_allele %in% alleles_tested, ras_allele, "Other"),
      allele = str_remove(allele, "KRAS_"),
      allele = factor_alleles(allele)
    )
  mod_stats_df <- stats_df %>%
    mutate(
      group1 = str_remove(group1, "KRAS_"),
      group2 = str_remove(group2, "KRAS_")
    )

  ridge_plots <- ggridge_plot_signatures_pgridge(
    mod_ms_df,
    signature = signature
  )

  stats_plot <- ggridge_plot_signatures_pstats(
    mod_stats_df,
    alleles = as.character(sort(unique(mod_ms_df$allele))),
    added_height = stats_plot_add_height
  )

  patch <- (stats_plot | ridge_plots) +
    plot_layout(widths = patch_widths)


  fn <- as.character(glue("ggrdige-stats_{cancer}_sig{signature}.svg"))
  ggsave_wrapper(
    patch,
    plot_path(GRAPHS_DIR, fn),
    "wide"
  )
}


plot_vars <- tribble(
  ~signature, ~cancer, ~stats_plot_add_height, ~patch_widths,
  "8", "COAD", 1.5, c(1, 10),
  "17", "COAD", 1.5, c(1, 10),
  "18", "COAD", 1.5, c(1, 10),
  "N3V2", "COAD", 1.5, c(1, 10),
  "2", "LUAD", 1.5, c(1, 10),
  "4", "LUAD", 1.5, c(1, 6),
  "15", "LUAD", 1.5, c(1, 10),
  "18", "LUAD", 1.5, c(1, 10),
  "2", "MM", 1.5, c(1, 10),
  "9", "MM", 1.5, c(1, 10),
  "1", "PAAD", 1.5, c(1, 10),
  "5", "PAAD", 1.5, c(1, 10),
  "6", "PAAD", 1.5, c(1, 10),
  "8", "PAAD", 1.5, c(1, 10),
  "17", "PAAD", 1.5, c(1, 10),
)


allele_signature_associations %>%
  filter(p_value < 0.05) %>%
  distinct(signature, cancer, group1, group2, p_value) %>%
  group_by(cancer, signature) %>%
  nest() %>%
  ungroup() %>%
  left_join(alleles_tested_per_cancer, by = "cancer") %>%
  rename(stats_df = data) %>%
  left_join(plot_vars, by = c("cancer", "signature")) %>%
  pwalk(ggridge_plot_signatures, ms_df = mutsig_noartifact_df)
