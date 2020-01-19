# A ggplot box-plot with comparison bars.

stats_boxplot <- function(df, stats_df,
                          box_color = NULL, box_fill = NULL,
                          box_alpha = 0.6, point_alpha = 0.9) {
    check_data(df)
    check_stats(stats_df)

    bp <- ggplot(df, aes(x = x, y = y)) +
        geom_boxplot(
            aes(color = box_color, fill = box_fill),
            alpha = box_alpha
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

    return(bp)
}


GRAPHS_DIR <- "TEST_STATS_BOXPLOT"
reset_graph_directory(GRAPHS_DIR)

set.seed(0)
test_df <- tibble(
    x = rep(LETTERS[1:3], 5),
    y = rnorm(length(x)),
    box_color = x
)

test_stats_df <- tibble(
    x1 = "A",
    x2 = "B", 
    label = "**"
)

ggsave_wrapper(
    stats_boxplot(test_df, test_stats_df),
    plot_path(GRAPHS_DIR, "test_box.svg"),
    "small"
)

check_data <- function(df) {
    assertr::verify(df, has_all_names("x", "y"))
}

check_stats <- function(df) {
    assertr::verify(df, has_all_names("x1", "x2", "label"))
}
