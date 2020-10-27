
# A ggplot box-plot with comparison bars.
#
# Create a box-plot with statistical comparison bars placed at the optimal
# location on the plot. They will not overlap with the data or each other. The
# function returns a 'ggplot2' box-plot without any custom components. The
# box-plot is made with `geom_boxplot(aes(x = x, y = y))` and the bars and
# labels are `annotate("text", ...)` objects.
#
# df: A data frame with the columns `x`, `y` to be plotted using
#   `geom_boxplot()`.
# stats_df: A data frame with the results of a statistical analysis between the
#   groups on the x-axis. It must have the column names `x1`, `x2`, and `label`
#   where the `label` is the result of the comparison of `x1` and `x2`. A bar
#   is drawn for each row of the data frame, so only include the comparisons
#   that should be on the final plot. If left `NULL`, one will be made using
#   the function `make_stats_dataframe()` using `auto_filter = TRUE`,
#   `method = "t.test"`, and `p.adjust.method = "BH"`.
# box_color, box_fill: passed to `color` and `fill` for `geom_boxplot(aes())`
#   and `geom_jitter(aes())`. If the properties are to be paired to a column
#   of the data frame, pass them unquoted (like a normal call to 'ggplot' for
#   setting an aesthetic). If a specific color/fill is to be used, pass it as a
#   string (e.g. "green") and apply `scale_color/fill_identity()`.
# point_alpha, point_shape, point_size, jitter_width, jitter_height: passed to
#  `alpha`, `shape`, `size`, `width`, and `height` for `geom_jitter()`.
# up_spacing, dn_spacing: The spacing above and below each bar. The units are
#   the same as the y-axis. Tweak these to provide enough space between the
#   bars for the labels.
# bar_color, bar_alpha, bar_size: The `color`, `alpha`, and `size` for
#   `geom_segment`.
# label_color, label_alpha, label_size, fontface, family: The `color`, `alpha`,
#   `size`, `fontface`, and `family` passed to `annotate("text")` for the
#   labels.
stats_boxplot <- function(df, stats_df = NULL,
                          box_color = NULL, box_fill = NULL,
                          box_alpha = 0.6, point_alpha = 0.9,
                          point_shape = NULL, point_size = 1,
                          jitter_width = NULL, jitter_height = NULL,
                          up_spacing = 0.1, dn_spacing = up_spacing,
                          bar_color = "black", bar_alpha = 1.0,
                          bar_size = 0.7, label_color = "black",
                          label_alpha = 1.0,
                          label_size = 8, fontface = "plain",
                          family = "Arial") {
  check_data(df)

  if (!is.null(stats_df)) {
    check_stats(stats_df)
  } else {
    stats_df <- make_stats_dataframe(df,
      auto_filter = TRUE,
      method = "t.test",
      p.adjust.method = "BH"
    )
  }

  bp <- ggplot(df, aes(x = x, y = y)) +
    geom_boxplot(
      aes(
        color = !!rlang::enquo(box_color),
        fill = !!rlang::enquo(box_fill)
      ),
      alpha = box_alpha, outlier.shape = NA
    ) +
    geom_jitter(
      aes(
        color = !!rlang::enquo(box_color),
        fill = !!rlang::enquo(box_fill)
      ),
      alpha = point_alpha, size = point_size,
      width = jitter_width, height = jitter_height
    )

  bp <- add_stats_comparisons(bp, df, stats_df,
    up_spacing = up_spacing,
    dn_spacing = dn_spacing,
    bar_color = bar_color,
    bar_alpha = bar_alpha,
    bar_size = bar_size,
    label_color = label_color,
    label_alpha = label_alpha,
    label_size = label_size,
    fontface = fontface,
    family = family
  )
  return(bp)
}


# Make a data frame that complies with the `stats_df` for `stats_boxplot()`.
#
# It must comply to the standards of that argument: have an `x` and `y` col.
# The `...` are passed to `ggpubr::compare_means()`.
#
# df: The data frame to be plotted. It must comply to the same standards as
#    described in `stats_boxplot()`.
# auto_filter: Should the results of the statistical analysis be filtered by
#   `label != "ns"`?
# ...: Passed to `ggpubr::compare_means()`.
make_stats_dataframe <- function(df, auto_filter = FALSE, ...) {
  check_data(df)

  stats_df <- compare_means(y ~ x, data = df, ...) %>%
    dplyr::rename(x1 = group1, x2 = group2, label = p.signif)

  if (auto_filter) {
    stats_df %<>% filter(label != "ns")
  }

  return(stats_df)
}


################################################################################
# GRAPHS_DIR <- "TEST_STATS_BOXPLOT"
# reset_graph_directory(GRAPHS_DIR)

# for (i in seq(1, 5)) {

#     set.seed(i)
#     n_test_samples <- 10
#     n_test_groups <- 5
#     test_df <- tibble(
#         x = unlist(purrr::map(1:n_test_groups,
#                               ~rep(LETTERS[.x], n_test_samples))),
#         y = unlist(purrr::map(runif(n_test_groups, 0, 5),
#                               ~rnorm(n_test_samples, mean = .x))),
#         box_color = x
#     )



#     ggsave_wrapper(
#         stats_boxplot(test_df, up_spacing = 0.3),
#         plot_path(GRAPHS_DIR, glue("test_box_{i}.svg")),
#         "small"
#     )
# }
################################################################################


# Check that the data frame `df` has the necessary columns.
check_data <- function(df) {
  assertr::verify(df, has_all_names("x", "y"))
}


# Check that the statistics data frame has the necessary columns.
check_stats <- function(df) {
  assertr::verify(df, has_all_names("x1", "x2", "label"))
}


# Create a new bar object.
make_new_bar <- function(x1, x2, y, up_space = 0, dn_space = up_space) {
  return(list(
    x1 = x1,
    x2 = x2,
    y = y,
    up_space = up_space,
    dn_space = dn_space,
    y_up = y + up_space,
    y_dn = y - dn_space
  ))
}


# Updeate the range of the bar.
update_bar <- function(bar) {
  bar$y_up <- bar$y + bar$up_space
  bar$y_dn <- bar$y - bar$dn_space
  return(bar)
}


# Is `x` between `left` and `right? Set inclusive to use greater/less than *or
# equal to*. `x` can be a vector, but `left` and `right` must be atomic values.
is_between <- function(x, left, right, inclusive = FALSE) {
  between_non_inclusive <- function(i) {
    return(left < i & i < right)
  }

  between_inclusive <- function(i) {
    return(left <= i & i <= right)
  }

  if (inclusive) {
    return(purrr::map_lgl(x, between_inclusive))
  } else {
    return(purrr::map_lgl(x, between_non_inclusive))
  }
}


# Do two bar objects overlap?
bars_overlap <- function(b1, b2) {
  # Check x-values first.
  x_cond1 <- is_between(b1$x1, b2$x1, b2$x2, inclusive = TRUE)
  x_cond2 <- is_between(b1$x2, b2$x1, b2$x2, inclusive = TRUE)
  if (!x_cond1 & !x_cond2) {
    return(FALSE)
  }

  # Check y-values if x-values overlap
  y_cond1 <- is_between(b1$y_dn, b2$y_dn, b2$y_up)
  y_cond2 <- is_between(b1$y_up, b2$y_dn, b2$y_up)
  y_cond3 <- is_between(b1$y, b2$y_dn, b2$y_up)
  return(any(c(y_cond1, y_cond2, y_cond3)))
}


# Move `b1` above `b2`.
move_b1_above_b2 <- function(b1, b2, buffer = 1e-4) {
  b1$y <- b2$y + b2$up_space + b1$dn_space + buffer
  b1 <- update_bar(b1)
  return(b1)
}


# Add the statistical comparison bars in `stats_df$bar` to the 'ggplot2'
# box-plot `bp`.
# See `stats_boxplot()` for information on the arguments.
add_stats_comparisons <- function(bp, df, stats_df,
                                  up_spacing = 0.1, dn_spacing = up_spacing,
                                  bar_color = "black", bar_alpha = 1.0,
                                  bar_size = 0.7, label_color = "black",
                                  label_alpha = 1.0,
                                  label_size = 8, fontface = "plain",
                                  family = "Arial") {
  stats_df %<>% arrange(x1, x2, label) %>%
    add_column(bar = NA) %>%
    mutate(
      x1_num = reassign_x(x1, list(x1, x2)),
      x2_num = reassign_x(x2, list(x1, x2))
    )

  df$x_num <- reassign_x(df$x, df$x)

  stats_df <- get_bars(df, stats_df, up_spacing, dn_spacing)

  for (i in seq(1, nrow(stats_df))) {
    bp <- add_bar_to_plot(bp, stats_df,
      row = i,
      color = bar_color, alpha = bar_alpha,
      size = bar_size
    )
    bp <- add_labels_to_plot(bp, stats_df,
      row = i,
      color = label_color, alpha = label_alpha,
      size = label_size, fontface = fontface,
      family = family
    )
  }

  max_y_val <- max(purrr::map_dbl(stats_df$bar, ~ unlist(.x)["y_up"]))
  bp <- bp +
    scale_y_continuous(limits = c(NA, max_y_val))

  return(bp)
}


# Assign `x` to a numeric value in the context of `all_xs`. If `x` is a factor,
# then the level is returned. Else, `x` are turned into a factor with
# levels determined by `sort(unique(x_vals))` and the level of `x` is returned.
# The `...` is passed to `sort()`.
reassign_x <- function(x, all_xs, ...) {
  if (is.factor(x)) {
    x <- forcats::fct_drop(x)
    return(as.integer(x))
  }
  if (!is.factor(x)) {
    fx <- factor(x, levels = sort(unique(unlist(all_xs)), ...))
    return(as.integer(fx))
  }
}


# A bar object is added to a plot. The bar should be in row `row` of data frame
# `stats_df`. The arguments `color`, `alpha`, and `size` are passed to
# `ggplot2::geom_segment()`.
add_bar_to_plot <- function(p, stats_df, row,
                            color = "black", alpha = 1.0, size = 0.7) {
  from <- stats_df$x1_num[[row]]
  to <- stats_df$x2_num[[row]]
  bar <- unlist(stats_df$bar[[row]], recursive = FALSE)
  p <- p +
    geom_segment(
      x = from, xend = to,
      y = bar$y, yend = bar$y,
      color = color, alpha = alpha, size = size
    )
  return(p)
}


# Statistical significance labels are added to the plot. The bar should be in
# row `row` of data frame `stats_df`. The arguments `color`, `alpha`, `size`,
# `fontface`, and `family` are passed to `ggplot2::annotate("text")`.
add_labels_to_plot <- function(p, stats_df, row,
                               color = "black", alpha = 1.0, size = 8,
                               fontface = "plain", family = "Arial") {
  from <- stats_df$x1_num[[row]]
  to <- stats_df$x2_num[[row]]
  mid <- mean(c(from, to))
  bar <- unlist(stats_df$bar[[row]], recursive = FALSE)

  lbl <- stats_df$label[[row]]
  if (is.na(lbl)) {
    return(p)
  }

  p <- p +
    annotate("text",
      x = mid, y = bar$y, label = lbl,
      colour = color, alpha = alpha, size = size,
      fontface = fontface, family = family
    )
  return(p)
}


# Find the optimal location for the bars such that they do not overlap with
# data nor each other. The spacing can be adjusted so that the labels do not
# overlap. The algorithm is as follows:
#  1. Place the first bar as low as it can go without overlapping any data.
#  2. Place the next bar as low as it can go without overlapping any data.
#  3. Check if the new bar overlaps any previously placed bars.
#  4. If it does, move it to have its lower bound just above the upper bound
#     of the lowest bar it overlaps with.
#  5. Repeat steps 3 and 4 until the new bar does not overlap with any others.
#  6. Repeat steps 2-5 until all bars have been placed.
get_bars <- function(df, stats_df, up_spacing, dn_spacing) {
  for (i in seq(1, nrow(stats_df))) {
    y_min <- lowest_y(stats_df$x1_num[[i]], stats_df$x2_num[[i]], df)
    b <- make_new_bar(
      x1 = stats_df$x1_num[[i]],
      x2 = stats_df$x2_num[[i]],
      y = y_min + dn_spacing,
      up_space = up_spacing,
      dn_space = dn_spacing
    )
    current_bars <- unlist(stats_df$bar, recursive = FALSE)
    b <- find_optimal_bar_position(b, current_bars)
    stats_df$bar[[i]] <- list(b)
  }
  return(stats_df)
}


# Find the lowest y-value for a bar.
lowest_y <- function(x1, x2, df) {
  xmin <- min(c(x1, x2))
  xmax <- max(c(x1, x2))

  min_vals <- df %>%
    filter(is_between(x_num, !!xmin, !!xmax, inclusive = TRUE)) %>%
    pull(y) %>%
    max()
  return(min_vals[[1]])
}


# Find the lowest bar in a list of bars.
find_lowest_bar <- function(bars) {
  ys <- purrr::map_dbl(bars, ~ .x$y)
  return(bars[[which.min(ys)]])
}


# Find the optimial location for a bar `b` amongst a list of bars `current_bs`.
# See `get_bars()` for the complete algorithm used.
find_optimal_bar_position <- function(b, current_bs) {
  current_bs <- current_bs[!is.na(current_bs)]

  if (length(current_bs) == 0) {
    return(b)
  }

  conflicts <- rep(TRUE, length(current_bs))

  while (any(conflicts)) {
    for (j in seq(1, length(conflicts))) {
      b2 <- current_bs[[j]]
      conflicts[[j]] <- bars_overlap(b, b2)
    }

    if (!any(conflicts)) break

    lowest_other_b <- find_lowest_bar(current_bs[conflicts])
    b <- move_b1_above_b2(b, lowest_other_b)
  }

  return(b)
}
