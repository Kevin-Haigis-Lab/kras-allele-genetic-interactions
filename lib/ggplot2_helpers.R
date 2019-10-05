
#' Get the full path for a plot name.
plot_path <- function(...) {
    file.path("graphs", ...)
}

#' Predetermined plot dimensions.
plot_size_dict <- list(
    small = c(4, 4),
    medium = c(6, 6),
    large = c(8, 8),
    tall = c(4, 8),
    wide = c(8, 4)
)

#' Save `ggplot` objects in a standardized fashion.
ggsave_wrapper <- function(p, save_path, size = NA, width = NA, height = NA) {
    size <- size[[1]]
    if (size %in% names(plot_size_dict)) {
        dims <- unlist(plot_size_dict[size])
    } else if (!is.na(width) & !is.na(height)) {
        dims <- c(width, height)
    } else if (!is.na(width)) {
        dims <- c(width, width)
    } else if (!is.na(height)) {
        dims <- c(height, height)
    } else {
        stop("Either select or input a width and height (or just one for a square)")
    }

    if (tools::file_ext(save_path) != "svg") {
        message("File name does not have a 'svg' extension.")
    }

    ggsave(filename = save_path, plot = p,
           width = dims[[1]], height = dims[[2]])
}
