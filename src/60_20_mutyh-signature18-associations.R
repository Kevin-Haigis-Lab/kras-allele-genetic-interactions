# Specifcially inspect the status of MUTYH in tumors with high Signature 18.
# This analysis was restricted to COAD and PAAD were Sig 18 had a relatively
# high probability of causing G12C mutations.

GRAPHS_DIR <- "60_20_mutyh-signature18-associations"
reset_graph_directory(GRAPHS_DIR)



#### ---- Prepare Signature 18 data ---- ####

sig18_df <- mutsig_noartifact_df %>%
  filter(signature == "18" & cancer %in% c("COAD", "PAAD"))


#### ---- Plot levels of Sig. 18 in samples ---- ####

plot_signature_density <- function(df) {
  p <- df %>%
    ggplot(aes(x = contribution)) +
    facet_wrap(~cancer, ncol = 1, scales = "free") +
    geom_density(aes(color = cancer, fill = cancer), alpha = 0.7) +
    scale_fill_manual(values = cancer_palette, guide = FALSE) +
    scale_color_manual(values = cancer_palette, guide = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      strip.background = element_blank()
    ) +
    labs(
      x = "signature 18 levels"
    )
  invisible(p)
}

sig18_density_plot_all <- sig18_df %>%
  plot_signature_density() +
  labs(title = "All samples")
sig18_density_plot_nozeros <- sig18_df %>%
  filter(contribution > 0.0) %>%
  plot_signature_density() +
  labs(title = "Without zeros")
sig18_density_plot_g12c <- sig18_df %>%
  filter(ras_allele == "KRAS_G12C") %>%
  plot_signature_density() +
  labs(title = "Without zeros")

sig18_density_plot <- (
  sig18_density_plot_all |
    sig18_density_plot_nozeros |
    sig18_density_plot_g12c
)

ggsave_wrapper(
  sig18_density_plot,
  plot_path(GRAPHS_DIR, "sig18_density_plot.svg"),
  width = 12, height = 4
)



plot_signature_waterfall <- function(df) {
  p <- df %>%
    mutate(tumor_sample_barcode = fct_reorder(
      tumor_sample_barcode,
      contribution
    )) %>%
    ggplot(aes(x = tumor_sample_barcode, y = contribution)) +
    facet_wrap(~cancer, scales = "free", ncol = 1) +
    geom_point(aes(color = colors, alpha = alphas), size = 0.5) +
    scale_color_identity() +
    scale_alpha_identity() +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(x = "ranked tumor samples", y = "sig. 18 level")
  invisible(p)
}


sig18_waterfall_plot <- sig18_df %>%
  mutate(
    colors = ifelse(ras_allele == "KRAS_G12C", "blue", "grey70"),
    alphas = ifelse(ras_allele == "KRAS_G12C", 1.0, 0.5)
  ) %>%
  plot_signature_waterfall()

ggsave_wrapper(
  sig18_waterfall_plot,
  plot_path(GRAPHS_DIR, "sig18_waterfall_plot.svg"),
  width = 12, height = 4
)


#### ---- Cluster using density-based method ---- ####


# Calculate the density of YFP values (log and scaled).
measure_density <- function(x,
                            log_scale = FALSE,
                            low_quantile = 0.001,
                            upper_quantile = 0.999,
                            n_bins = 500,
                            min_density = 0.0) {
  if (log_scale) {
    x <- scale(log(x))
  }

  # Remove outliers by filtering on quantiles.
  lower_cut <- quantile(x, low_quantile, na.rm = TRUE)
  upper_cut <- quantile(x, upper_quantile, na.rm = TRUE)
  x <- x[x > lower_cut & x < upper_cut]

  d <- density(x, n = n_bins, na.rm = TRUE)


  # Only include intervals above a minimum density
  idx <- (d$y > min_density)

  return(list(dx = d$x[idx], dy = d$y[idx]))
}


# Is the center value in `vals` the maximum value?
center_is_maximum <- function(vals) {
  mid_point <- ceiling(length(vals) / 2)
  return(which.max(vals) == mid_point)
}


# Locate local maxima.
find_local_max <- function(dx, dy, n_either_side = 3) {
  start_i <- n_either_side
  finish_i <- length(dx) - n_either_side
  local_maxes <- c()
  for (i in seq(start_i, finish_i)) {
    vals <- dy[seq((i - n_either_side), (i + n_either_side))]
    if (center_is_maximum(vals)) {
      local_maxes <- c(local_maxes, dx[[i]])
    }
  }
  return(sort(local_maxes))
}


# Locate local minima between the local maxima.
# If `maxes` if `NULL`, the local maxima are calculated.
# There must be at least 2 local maxima to find local minima.
find_local_min <- function(x, y, maxes = NULL) {
  if (is.null(maxes)) maxes <- find_local_max(x, y)

  if (length(maxes) < 2) {
    warning("Not enough local maxima to find local minima.")
    return(NULL)
  }

  local_mins <- c()
  for (i in seq(1, length(maxes) - 1)) {
    m_low <- maxes[[i]]
    m_high <- maxes[[i + 1]]

    idx <- (x > m_low & x < m_high)
    x_interval <- x[idx]
    y_interval <- y[idx]

    interval_min <- which.min(y_interval)
    local_mins <- c(local_mins, x_interval[interval_min])
  }

  return(sort(local_mins))
}


add_zero_padding <- function(df, num, min_val = -0.03, max_val = 0) {
  new_rows <- head(df, num)
  new_rows$contribution <- seq(min_val, max_val, length.out = num)
  return(bind_rows(df, new_rows))
}


classify_by_density <- function(df,
                                zero_padding = NULL,
                                slide_width = 5,
                                ...) {
  set.seed(0)

  if (!is.null(zero_padding)) {
    df <- add_zero_padding(df, num = zero_padding)
  }

  densities <- measure_density(df$contribution, ...)
  local_maxes <- find_local_max(densities$dx, densities$dy,
    n_either_side = slide_width
  )
  local_mins <- find_local_min(densities$dx, densities$dy,
    maxes = local_maxes
  )

  mdl <- list(
    densities = densities,
    local_maxima = local_maxes,
    local_minima = local_mins
  )

  if (length(local_mins) == 0) {
    df$classification <- 0
  }

  df$classification <- length(local_mins) + 1
  for (i in length(local_mins)) {
    df$classification[df$contribution < local_mins[[i]]] <- i
  }

  return(list(data = df, model = mdl))
}


plot_classified_signature_density <- function(data, ...) {
  p <- data %>%
    ggplot(aes(x = contribution)) +
    facet_wrap(~cancer, ncol = 1, scales = "free") +
    geom_density(
      aes(color = classification, fill = classification),
      alpha = 0.7
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      strip.background = element_blank()
    ) +
    labs(
      x = "signature 18 levels"
    )
  invisible(p)
}


plot_density_model <- function(density_model, cancer, ...) {
  mdl <- density_model
  p <- tibble(x = mdl$densities$dx, y = mdl$densities$dy) %>%
    mutate(pts = case_when(
      x %in% mdl$local_maxima ~ "maxima",
      x %in% mdl$local_minima ~ "minima"
    )) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(color = "black") +
    geom_point(aes(color = pts), size = 0.7) +
    scale_color_manual(values = c("dodgerblue", "tomato")) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "x-values", y = "density", title = cancer)
  invisible(p)
}



sig18_df_cls <- sig18_df %>%
  filter(contribution > 0.0) %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(
    density_classes = purrr::map(data, classify_by_density,
      n_bins = 500,
      low_quantile = 0.0,
      upper_quantile = 1.0
    ),
    data = purrr::map(density_classes, ~ .x$data),
    density_model = purrr::map(density_classes, ~ .x$model)
  ) %>%
  select(-density_classes)

classified_sig18_density_plot <- sig18_df_cls %>%
  select(cancer, data) %>%
  unnest(data) %>%
  mutate(classification = as.character(classification)) %>%
  plot_classified_signature_density()

ggsave_wrapper(
  classified_sig18_density_plot,
  plot_path(GRAPHS_DIR, "classified_sig18_density_plot.svg"),
  "small"
)

model_density_plots <- sig18_df_cls %>%
  purrr::pmap(plot_density_model)
model_density_plots <- wrap_plots(model_density_plots, ncol = 1)
ggsave_wrapper(
  model_density_plots,
  plot_path(GRAPHS_DIR, "model_density_plots.svg"),
  "small"
)



sig18_highs <- sig18_df %>%
  filter(contribution > 0.0) %>%
  select(cancer, tumor_sample_barcode, ras_allele) %>%
  unique()

table(sig18_highs$cancer)
table(sig18_highs$cancer, sig18_highs$ras_allele)


#### ---- MUTYH status of Sig. 18-high G12C tumors ---- ####

cancer_full_muts_df %>%
  filter(tumor_sample_barcode %in% sig18_highs$tumor_sample_barcode) %>%
  filter(hugo_symbol == "MUTYH")
