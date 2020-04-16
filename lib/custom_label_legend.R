# A custom lagend using `geom_label()`.
#
# limitation: only works for a single row of labels
#
# Use the separate functions to manually tweak any values.
custom_label_legend <- function(lbl,
                                gap = 0,
                                colors = rep("black", length(lbl)),
                                ...) {
    custom_label_legend_df(lbl, gap, colors) %>%
        custom_label_legend_plot(...)
}

custom_label_legend_df <- function(lbl,
                                   gap = 0,
                                   colors = rep("black", length(lbl))) {
    tibble(lbl = lbl) %>%
        mutate(
            len = get_string_length(lbl),
            start = calc_starts(len, gap = gap),
            end = len + start,
            mid = (start + end) / 2
        ) %>%
        add_column(color = colors)
}


custom_label_legend_plot <- function(legend_df, ...) {
    legend_df %>%
        ggplot(aes(x = mid, y = 1, label = lbl, fill = lbl, color = color)) +
        geom_label(...) +
        scale_color_identity() +
        scale_x_continuous(
            limits = c(min(legend_df$start), max(legend_df$end)),
            expand = c(0, 0)
        ) +
        theme_void()
}


calc_starts <- function(x, gap = 0, starting_pt = 0) {
    starts <- accumulate(x, .init = starting_pt, ~ .x + .y + gap)
    starts[-length(starts)]
}


get_string_length <- function(x) {
    if (all(x %in% string_length_dictionary$word)) {
        known_lengths <- tibble(word = x) %>%
            left_join(string_length_dictionary, by = "word") %>%
            pull(known_length) %>%
            unlist()
        return(known_lengths / min(known_lengths))
    } else {
        return(str_length(x))
    }
}


string_length_dictionary <- tibble::tribble(
    ~word, ~known_length,
    "missense", 0.7039388,
    "frame shift del.", 1.0930176,
    "frame shift ins.", 1.0836589,
    "nonsense", 0.7228190,
    "in-frame ins.", 0.9169108,
    "splice site", 0.7317708
)



calc_string_width <- function(s) {
    strwidth(s, font = 12, units = "in", family = "Arial")
}