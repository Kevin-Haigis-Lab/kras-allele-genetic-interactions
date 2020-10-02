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


theme_set(
  theme_bw(base_size = 11, base_family = "Arial") %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.ticks = element_blank()
    )
)


allele_signature_associations <- readRDS("~/Downloads/allele_signature_associations.RDS")
allele_signature_causation_stats <- readRDS("~/Downloads/allele_signature_causation_stats.RDS")
MOD_kras_allele_causation_mutsig_df <- readRDS("~/Downloads/MOD_kras_allele_causation_mutsig_df.rds")
mutsig_noartifact_df <- readRDS("~/Downloads/mutsig_noartifact_df.rds")


################################################################################
# COPPIED

GRAPHS_DIR <- "50_35_mutational-signature-allele-associations"


modify_mutsig_names <- function(sig) {
  case_when(
    sig == "N3V2" ~ "N",
    TRUE ~ sig
  )
}


factor_signatures <- function(sig) {
  factor(sig, levels = c(names(mutsig_pal), "N"))
}


strip_ras <- function(df, col = ras_allele) {
  df %>%
    mutate({{ col }} := str_remove({{ col }}, "KRAS_"))
}


prep_mutsig_dataframe_for_plotting <- function(ms_df,
                                               cancer,
                                               signature,
                                               alleles_tested) {
  ms_df %>%
    strip_ras(ras_allele) %>%
    mutate(signature = modify_mutsig_names(signature)) %>%
    filter(cancer == !!cancer & signature == !!signature) %>%
    mutate(
      allele = ifelse(ras_allele %in% !!alleles_tested, ras_allele, "Other"),
      allele = factor_alleles(allele),
      allele = fct_drop(allele)
    )
}



################################################################################




boxplot_plot_signatures_pboxes <- function(df, signature, y_title = NULL) {
  if (is.null(y_title)) {
    y_title <- "mutational signature {signature}"
  }

  df %>%
    ggplot(aes(x = allele, y = contribution)) +
    geom_boxplot(
      color = "grey25",
      fill = "grey50",
      alpha = 0.5,
      width = 0.5,
      outlier.shape = NA
    ) +
    geom_jitter(
      alpha = 0.3,
      size = 0.3,
      width = 0.2,
      height = 0
    ) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      panel.grid.major.x = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      x = NULL,
      y = glue(y_title),
      title = NULL
    )
}


get_max_value_between <- function(g1, g2, df) {
  n1 <- as.numeric(g1)
  n2 <- as.numeric(g2)
  df %>%
    filter(between(as.numeric(allele), min(n1, n2), max(n1, n2))) %>%
    pull(contribution) %>%
    unlist() %>%
    max()
}


optimize_box_stats_bar_placement <- function(df, raise_by) {
  for (i in sort(df$idx)) {
    if (i == 1) {
      next
    }

    g1_i <- as.numeric(df$group1[[i]])
    g2_i <- as.numeric(df$group2[[i]])
    n1_i <- min(g1_i, g2_i)
    n2_i <- max(g1_i, g2_i)

    is_overlaping <- TRUE
    k <- 0

    while (is_overlaping) {
      k <- k + 1
      if (k > 20) {
        stop("Loop has continued for 10 iterations!")
      }

      y_i <- df$y[[i]]

      other_df <- df %>%
        filter(idx < i & near(y, y_i, tol = 0.01))

      if (nrow(other_df) == 0) {
        is_overlaping <- FALSE
        break
      }

      other_df <- other_df %>%
        filter(
          between(as.numeric(group1), n1_i, n2_i) |
            between(as.numeric(group2), n1_i, n2_i) |
            (
              as.numeric(group1) <= n1_i & as.numeric(group1) <= n2_i &
                as.numeric(group2) >= n1_i & as.numeric(group2) >= n2_i
            ) |
            (
              as.numeric(group2) <= n1_i & as.numeric(group2) <= n2_i &
                as.numeric(group1) >= n1_i & as.numeric(group1) >= n2_i
            )
        )

      if (nrow(other_df) == 0) {
        is_overlaping <- FALSE
        break
      }

      df$y[i] <- df$y[i] + raise_by
    }
  }
  return(df)
}


create_boxplot_stats_plotting_dataframe <- function(ms_df, stats_df, fct_dist = 0.05) {
  raise_by <- fct_dist * (max(ms_df$contribution) - min(ms_df$contribution))

  calc_sort_value <- function(g1, g2) {
    g1 <- as.numeric(g1)
    g2 <- as.numeric(g2)
    n1 <- min(g1, g2)
    n2 <- max(g1, g2)
    return(n1 * n2 * n2)
  }

  stats_df %>%
    mutate(sort_value = map2_dbl(group1, group2, calc_sort_value)) %>%
    arrange(sort_value) %>%
    mutate(
      idx = row_number(),
      y = map2_dbl(group1, group2, get_max_value_between, df = ms_df),
      y = y + raise_by
    ) %T>%
    print() %>%
    optimize_box_stats_bar_placement(raise_by = raise_by) %T>%
    print()
}


add_star_label_data <- function(stats_df) {
  stats_df %>%
    mutate(
      stars_lbl = assign_stars(p_value),
      stars_x = map2_dbl(
        as.numeric(group1),
        as.numeric(group2),
        ~ mean(c(.x, .y))
      )
    )
}


annotate_boxplot_with_statbars <- function(bp, bars_df) {
  long_bars_df <- bars_df %>%
    select(idx, group1, group2, y) %>%
    pivot_longer(-c(idx, y), names_to = NULL, values_to = "x")

  bp +
    geom_line(
      aes(x = x, y = y, group = idx),
      data = long_bars_df,
      size = 0.5,
      color = "grey15"
    ) +
    geom_text(
      aes(x = stars_x, y = y, label = stars_lbl),
      vjust = 0.3,
      size = 3.5,
      family = "Arial",
      color = "grey15",
      data = bars_df
    )
}

mutsig_boxplot_breaks <- function(ms_df) {
  min_y <- 0
  max_y <- round(max(ms_df$contribution), 1)
  pretty(c(min_y, max_y), 4)
}


boxplot_plot_signatures <- function(signature,
                                    cancer,
                                    ms_df,
                                    stats_df,
                                    fn_glue,
                                    alleles_tested = NULL,
                                    box_y_title = NULL,
                                    y_expand_up = 0.04,
                                    ...) {
  mod_ms_df <- prep_mutsig_dataframe_for_plotting(
    ms_df,
    cancer = cancer,
    signature = signature,
    alleles_tested = alleles_tested
  )

  box_plot <- boxplot_plot_signatures_pboxes(
    mod_ms_df,
    signature = signature,
    y_title = box_y_title
  ) + scale_y_continuous(
    breaks = mutsig_boxplot_breaks(mod_ms_df),
    expand = expansion(mult = c(0, y_expand_up))
  )

  ras_levels <- as.character(unique(sort(mod_ms_df$allele)))

  mod_stats_df <- stats_df %>%
    mutate(
      group1 = factor(group1, levels = ras_levels),
      group2 = factor(group2, levels = ras_levels)
    )

  stats_bars_df <- create_boxplot_stats_plotting_dataframe(
    mod_ms_df,
    mod_stats_df
  ) %>%
    add_star_label_data()
  box_plot <- annotate_boxplot_with_statbars(box_plot, stats_bars_df)
  # ggsave_wrapper(
  #   box_plot,
  #   plot_path(GRAPHS_DIR, glue(fn_glue)),
  #   "small"
  # )
  return(box_plot)
}


# box_plots_levels <- allele_signature_associations %>%
#   pmap(
#     boxplot_plot_signatures,
#     ms_df = mutsig_noartifact_df,
#     fn_glue = "boxplot_sig-levels_stats_{cancer}_sig{signature}.svg",
#     box_y_title = "mutational signature {signature} composition"
#   )
# saveFigRds(box_plots_levels, "box_plots_sig_levels")

box_plots_cause <- allele_signature_causation_stats %>%
  slice(6) %>%
  pmap(
    boxplot_plot_signatures,
    ms_df = MOD_kras_allele_causation_mutsig_df,
    fn_glue = "boxplot_sig-causation_stats_{cancer}_sig{signature}.svg",
    box_y_title = "prob. of causation by signature {signature}"
  )
# saveFigRds(box_plots_cause, "box_plots_sig_cause")
