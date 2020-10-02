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
      outlier.shape = NA
    ) +
    geom_jitter(
      alpha = 0.4,
      size = 0.4,
      width = 0.3,
      height = 0
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.02))
    ) +
    theme_bw(base_size = 11, base_family = "Arial") +
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

    if (i == 1) { next }

    g1_i <- as.numeric(df$group1[[i]])
    g2_i <- as.numeric(df$group2[[i]])
    n1_i <- min(g1_i, g2_i)
    n2_i <- max(g1_i, g2_i)

    is_overlaping <- TRUE

    while(is_overlaping) {

      y_i <- df$y[[i]]

      other_df <- df %>%
        filter(idx < i & near(y, y_i))

      if(nrow(other_df) == 0) {
        is_overlaping <- FALSE
        break
      }

      other_df <- other_df %>%
        filter(
          between(as.numeric(group1), n1_i, n2_i) |
            between(as.numeric(group2), n1_i, n2_i)
        )

      if(nrow(other_df) == 0) {
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
  stats_df %>%
    mutate(sort_value = as.numeric(group1) * as.numeric(group2)) %>%
    arrange(sort_value) %>%
    mutate(
      y = map2_dbl(group1, group2, get_max_value_between, df = ms_df),
      y = y + raise_by
    ) %>%
    optimize_box_stats_bar_placement(raise_by = raise_by)
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



boxplot_plot_signatures <- function(signature,
                                    cancer,
                                    ms_df,
                                    stats_df,
                                    fn_glue,
                                    alleles_tested = NULL,
                                    stats_plot_add_height = 1.5,
                                    patch_heights = c(1, 10),
                                    box_y_title = NULL) {
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
  )

  ras_levels <- as.character(unique(sort(mod_ms_df$allele)))

  mod_stats_df <- stats_df %>%
    mutate(
      idx = row_number(),
      group1 = factor(group1, levels = ras_levels),
      group2 = factor(group2, levels = ras_levels)
    )

  stats_bars_df <- create_boxplot_stats_plotting_dataframe(mod_ms_df,
                                                           mod_stats_df) %>%
    add_star_label_data()
  annotate_boxplot_with_statbars(box_plot, stats_bars_df)
}

annotate_boxplot_with_statbars <- function(bp, bars_df) {
  long_bars_df <- bars_df %>%
    select(idx, group1, group2, y) %>%
    pivot_longer(-c(idx, y), names_to = NULL, values_to = "x")

  bp +
    geom_line(aes(x = x, y = y, group = idx), data = long_bars_df) +
    geom_text(aes(x = stars_x, y = y, label = stars_lbl), data = bars_df)
}



allele_signature_associations %>%
  rename(patch_heights = patch_widths) %>%
  slice(1) %>%
  pmap(
    boxplot_plot_signatures,
    ms_df = mutsig_noartifact_df,
    fn_glue = "boxplot_sig-levels_stats_{cancer}_sig{signature}.svg",
    box_y_title = "mutational signature {signature} composition"
  )

