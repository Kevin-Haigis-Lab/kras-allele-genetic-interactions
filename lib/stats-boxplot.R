# A ggplot box-plot with comparison bars.

stats_boxplot <- function(df, stats_df,
                          box_color = NULL, box_fill = NULL,
                          box_alpha = 0.6, point_alpha = 0.9,
                          point_shape = NULL, point_size = NULL,
                          jitter_width = NULL, jitter_height = NULL,
                          up_spacing = 0.1, dn_spacing = up_spacing,
                          bar_color = "black", bar_alpha = 1.0,
                          bar_size = 0.7, label_color = "black",
                          label_alpha = 1.0,
                          label_size = 8, fontface = "plain",
                          family = "Arial"
                      ) {
    check_data(df)
    check_stats(stats_df)

    bp <- ggplot(df, aes(x = x, y = y)) +
        geom_boxplot(
            aes(color = box_color, fill = box_fill),
            alpha = box_alpha, outlier.shape = NA
        ) +
        geom_jitter(
            aes(color = box_color, fill = box_fill),
            alpha = point_alpha
        )


    if (!is.null(box_color)) {
        bp <- bp + scale_color_identity()
    }

    if (!is.null(box_fill)) {
        bp <- bp + scale_fill_identity()
    }

    bp <- add_stats_comparisons(bp, df, stats_df,
                                up_spacing = 0.1, dn_spacing = up_spacing,
                                bar_color = bar_color, bar_alpha = bar_alpha,
                                bar_size = bar_size, label_color = label_color,
                                label_alpha = label_alpha,
                                label_size = label_size, fontface = fontface,
                                family = family)
    return(bp)
}


################################################################################
GRAPHS_DIR <- "TEST_STATS_BOXPLOT"
reset_graph_directory(GRAPHS_DIR)

for (i in seq(1, 30)) {

    set.seed(i)

    test_df <- tibble(
        x = rep(LETTERS[1:5], 6),
        y = rnorm(length(x)),
        box_color = x
    )

    test_stats_df <- combn(unique(test_df$x), 2) %>%
        t() %>%
        as_tibble() %>%
        dplyr::rename(x1 = V1, x2 = V2) %>%
        mutate(label = "**") %>%
        sample_frac(sample(seq(0.1, 1.0, 0.1), 1))

    ggsave_wrapper(
        stats_boxplot(test_df, test_stats_df),
        plot_path(GRAPHS_DIR, glue("test_box_{i}.jpeg")),
        "small"
    )
}

################################################################################




check_data <- function(df) {
    assertr::verify(df, has_all_names("x", "y"))
}


check_stats <- function(df) {
    assertr::verify(df, has_all_names("x1", "x2", "label"))
}


# Create a new bar object.
make_new_bar <- function(x1, x2, y, up_space = 0, dn_space = up_space) {
    return(list(x1 = x1,
                x2 = x2,
                y = y,
                up_space = up_space,
                dn_space = dn_space,
                y_up = y + up_space,
                y_dn = y - dn_space))
}


# Updeate the range of the bar.
update_bar <- function(bar) {
    bar$y_up <- bar$y + bar$up_space
    bar$y_dn <- bar$y - bar$dn_space
    return(bar)
}


is_between <- function(x, left, right) {
    f <- function(i) {
        return(left < i & i < right)
    }
    purrr::map_lgl(x, f)
}


# Do two bar objects overlap?
bars_overlap <- function(b1, b2) {
    # Check x-values first.
    x_cond1 <- between(b1$x1, b2$x1, b2$x2)
    x_cond2 <- between(b1$x2, b2$x1, b2$x2)
    if (!x_cond1 & !x_cond2) return(FALSE)

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


add_stats_comparisons <- function(bp, df, stats_df,
                                  up_spacing = 0.1, dn_spacing = up_spacing,
                                  bar_color = "black", bar_alpha = 1.0,
                                  bar_size = 0.7, label_color = "black",
                                  label_alpha = 1.0,
                                  label_size = 8, fontface = "plain",
                                  family = "Arial") {
    stats_df %<>% arrange(x1, x2, label) %>%
        add_column(bar = NA) %>%
        mutate(x1_num = reassign_x(x1, list(x1, x2)),
               x2_num = reassign_x(x2, list(x1, x2)))

    df$x_num <- reassign_x(df$x, df$x)

    stats_df <- get_bars(df, stats_df, up_spacing, dn_spacing)

    for (i in seq(1, nrow(stats_df))) {
        bp <- add_bar_to_plot(bp, stats_df, row = i,
                              color = bar_color, alpha = bar_alpha,
                              size = bar_size)
        bp <- add_labels_to_plot(bp, stats_df, row = i,
                                 color = label_color, alpha = label_alpha,
                                 size = label_size, fontface = fontface,
                                 family = family)
    }


    max_y_val <- max(purrr::map_dbl(stats_df$bar, ~ unlist(.x)["y_up"]))
    bp <- bp +
        scale_y_continuous(limits = c(NA, max_y_val))

    return(bp)
}


reassign_x <- function(x, all_xs) {
    if (is.factor(x)) {
        return(as.numeric(x))
    }
    if (!is.factor(x)) {
        fx <- factor(x, levels = sort(unique(unlist(all_xs))))
        return(as.numeric(fx))
    }
}


add_bar_to_plot <- function(p, stats_df, row,
                            color = "black", alpha = 1.0, size = 0.7) {
    from <- stats_df$x1_num[[row]]
    to <- stats_df$x2_num[[row]]
    bar <- unlist(stats_df$bar[[row]], recursive = FALSE)
    p <- p +
        geom_segment(x = from, xend = to,
                     y = bar$y, yend = bar$y,
                     color = color, alpha = alpha, size = size)
    return(p)
}


add_labels_to_plot <- function(p, stats_df, row,
                               color = "black", alpha = 1.0, size = 8,
                               fontface = "plain", family = "Arial") {
    from <- stats_df$x1_num[[row]]
    to <- stats_df$x2_num[[row]]
    mid <- mean(c(from, to))
    bar <- unlist(stats_df$bar[[row]], recursive = FALSE)

    lbl <- stats_df$label[[row]]
    if (is.na(lbl)) return(p)

    p <- p +
        annotate("text", x = mid, y = bar$y, label = lbl,
                 colour = color, alpha = alpha, size = size,
                 fontface = fontface, family = family)
    return(p)
}


get_bars <- function(df, stats_df, up_spacing, dn_spacing) {
    for (i in seq(1, nrow(stats_df))) {
        y_min <- lowest_y(stats_df$x1_num[[i]], stats_df$x2_num[[i]], df)
        b <- make_new_bar(x1 = stats_df$x1_num[[i]],
                          x2 = stats_df$x2_num[[i]],
                          y = y_min + dn_spacing,
                          up_space = up_spacing,
                          dn_space = dn_spacing)
        current_bars <- unlist(stats_df$bar, recursive = FALSE)
        b <- find_optimal_bar_position(b, current_bars)
        stats_df$bar[[i]] <- list(b)
    }
    return(stats_df)
}


lowest_y <- function(x1, x2, df) {
    min_vals <- df %>%
        filter(between(x_num, x1, x2)) %>%
        pull(y) %>%
        max()
    return(min_vals[[1]])
}


find_lowest_bar <- function(bars) {
    ys <- purrr::map_dbl(bars, ~.x$y)
    return(bars[[which.min(ys)]])
}


find_optimal_bar_position <- function(b, current_bs) {
    current_bs <- current_bs[!is.na(current_bs)]

    if (length(current_bs) == 0) return(b)

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
