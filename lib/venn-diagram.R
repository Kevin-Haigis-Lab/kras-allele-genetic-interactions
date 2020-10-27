library(mustashe)
library(jhcutils)
library(glue)
library(ggtext)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
  if (str_detect(f, "global|enrich|venn")) {
    next
  }
  source(f)
}

options(dplyr.summarise.inform = FALSE)

################################################################################
################################################################################




# Create ggplot2 based Venn diagrams.

library(ggforce)
library(magrittr)
library(tidyverse)

#### ---- Area of intersection ---- ####

# source:
#   https://diego.assencio.com/?index=8d6ca3d82151bad815f78addf9b5c1c6

# Calculate the area of the intersection of 2 circles with radii `r1` and `r2`
# that are `d` far away from each other.
area_intersection_circles <- function(d, r1, r2) {
  if (r1 < r2) {
    HOLD_r2 <- r2
    r2 <- r1
    r1 <- HOLD_r2
  }

  if (d >= r1 + r2) {
    return(0)
  }
  if (d <= r1 - r2) {
    return(pi * r2^2)
  }

  d1 <- (r1^2 - r2^2 + d^2) / (2 * d)
  d2 <- d - d1
  A <- r1^2 * acos(d1 / r1) - d1 * sqrt(r1^2 - d1^2)
  A <- A + r2^2 * acos(d2 / r2) - d2 * sqrt(r2^2 - d2^2)
  return(A)
}


# A plot to show plot the area of intersection over different values.
# (Currently not used, but useful for diagnostics.)
optimization_diagnositic_plot <- function(target_area, r1, r2, x, y) {
  pal <- c("red", "green", "blue")
  x_lbl <- c(r1, r2, abs(r1 - r2))

  min_area <- min(c(r1, r2))^2 * pi * 0.5

  tibble(x, y) %>%
    ggplot(aes(x, y)) +
    geom_vline(xintercept = x_lbl, color = pal) +
    geom_hline(yintercept = c(target_area, min_area), color = c("purple", "orange")) +
    annotate("text",
      x = x_lbl + 0.01, y = mean(y),
      label = c("r1", "r2", "r1 - r2"),
      color = pal,
      hjust = 0
    ) +
    geom_point(size = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "grey25")
    ) +
    labs(
      x = "distance between centers",
      y = "intersection area"
    )
}


# Find the distance between two circles with radii `r1` and `r2` such that
# their overlap has area `target_area`
distance_for_area_intersection <- function(target_area, r1, r2,
                                           tol = .Machine$double.eps^0.25) {
  f <- function(x) {
    abs(target_area - area_intersection_circles(x, r1 = r1, r2 = r2))
  }

  lower <- abs(r1 - r2) - 0.001
  upper <- r1 + r2 + 0.001

  # x <- seq(0, upper, 0.001)
  # y <- map_dbl(x, area_intersection_circles, r1 = r1, r2 = r2)
  # p <- optimization_diagnositic_plot(target_area, r1, r2, x, y)

  optimize(
    f = f,
    interval = c(lower, upper),
    tol = tol
  )$minimum
}


#### ---- Plotting functions ---- ####


radius_from_area <- function(A) {
  sqrt(A / pi)
}


make_circle_dataframe <- function(r1, r2, x2, circle_fill) {
  df <- tibble(
    x0 = c(0, x2),
    y0 = c(0, 0),
    r = c(r1, r2)
  )

  if (!is.null(circle_fill)) {
    df$fill <- circle_fill
  } else {
    df$fill <- "white"
  }

  return(df)
}


make_count_labels_dataframe <- function(r1, r2,
                                        x2,
                                        n1_diff, n2_diff,
                                        n_intersect) {
  count_labels_df <- tibble(
    x = c((-r1 + x2 - r2) / 2, (x2 - r2 + r1) / 2, (r1 + x2 + r2) / 2),
    y = 0,
    label = c(n1_diff, n_intersect, n2_diff)
  )
}


ggvenndiagram <- function(s1, s2,
                          cat_label = NULL,
                          cat_shift_up = 0.1,
                          cat_color = "black",
                          cat_size = NA,
                          cat_family = NA,
                          cat_fontface = "plain",
                          circle_fill = NULL,
                          circle_alpha = 0,
                          circle_color = "black",
                          circle_size = 0.6,
                          count_color = "black",
                          count_size = NA,
                          count_family = cat_family,
                          count_fontface = "plain",
                          count_labels_df = NULL) {
  n1 <- n_distinct(s1)
  n2 <- n_distinct(s2)
  n_intersect <- length(intersect(s1, s2))
  n1_diff <- length(setdiff(s1, s2))
  n2_diff <- length(setdiff(s2, s1))

  norm_fct <- length(union(s1, s2))
  a1 <- n1 / norm_fct
  a2 <- n2 / norm_fct
  a_intersect <- n_intersect / norm_fct

  r1 <- radius_from_area(a1)
  r2 <- radius_from_area(a2)

  x2 <- distance_for_area_intersection(
    target_area = a_intersect,
    r1 = r1,
    r2 = r2
  )

  circle_df <- make_circle_dataframe(r1, r2, x2, circle_fill)

  if (is.null(count_labels_df)) {
    count_labels_df <- make_count_labels_dataframe(
      r1, r2, x2, n1_diff, n2_diff, n_intersect
    )
  }

  count_labels_df %<>% filter(label > 0)

  p <- ggplot(circle_df) +
    geom_circle(
      aes(
        x0 = x0,
        y0 = y0,
        r = r,
        fill = fill
      ),
      color = circle_color,
      alpha = circle_alpha,
      size = circle_size
    ) +
    geom_text(
      aes(x = x, y = y, label = label),
      data = count_labels_df,
      color = count_color,
      size = count_size,
      family = count_family,
      fontface = count_fontface
    ) +
    scale_fill_identity()

  if (!is.null(cat_label)) {
    p <- p +
      annotate(
        "text",
        x = circle_df$x0,
        y = max(circle_df$r) + cat_shift_up,
        label = cat_label,
        color = cat_color,
        size = cat_size,
        family = cat_family,
        fontface = cat_fontface
      )
  }

  p <- p +
    labs(
      x = NULL,
      y = NULL
    )

  return(p)
}


#### ---- Examples ---- ####

if (FALSE) {
  venn_diagram_dir <- "TEST_VENN_DIAGRAM"
  reset_graph_directory(venn_diagram_dir)

  vd_1 <- ggvenndiagram(
    s1 = LETTERS[1:15],
    s2 = LETTERS[10:20],
    cat_label = c("group A", "group B")
  ) +
    coord_fixed() +
    theme_void()
  ggsave_wrapper(
    vd_1,
    plot_path(venn_diagram_dir, "venn_diagram_example_1.svg"),
    "small"
  )

  vd_2 <- ggvenndiagram(
    s1 = LETTERS[1:15],
    s2 = LETTERS[10:20],
    cat_label = c("group A", "group B"),
    cat_shift_up = -1,
    cat_color = "red",
    cat_size = 10,
    cat_family = "Times",
    cat_fontface = "bold",
    circle_fill = c("A", "B"),
    circle_alpha = 0.2,
    circle_color = "green",
    circle_size = 2,
    count_color = "purple",
    count_size = 20,
    count_family = "Arial",
    count_fontface = "italic",
    count_labels_df = NULL
  ) +
    scale_fill_discrete() +
    coord_fixed() +
    theme_void()
  ggsave_wrapper(
    vd_2,
    plot_path(venn_diagram_dir, "venn_diagram_example_2.svg"),
    "small"
  )

  vd_3 <- ggvenndiagram(
    s1 = LETTERS[1:15],
    s2 = LETTERS[16:26],
    cat_label = c("group A", "group B")
  ) +
    coord_fixed() +
    theme_void()
  ggsave_wrapper(
    vd_3,
    plot_path(venn_diagram_dir, "venn_diagram_example_3.svg"),
    "small"
  )

  # Remove unneeded variable.
  rm("venn_diagram_dir")
}
