
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


# decide on the sizes to use
decide_size <- function(size = NA, width = NA, height = NA) {
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
    return(dims)
}

#' Save `ggplot` objects in a standardized fashion.
ggsave_wrapper <- function(p, save_path, size = NA, width = NA, height = NA) {
    size <- decide_size(size = size[[1]], width = width, height = height)

    if (tools::file_ext(save_path) != "svg") {
        message("File name does not have a 'svg' extension.")
    }

    ggsave(filename = save_path, plot = p,
           width = dims[[1]], height = dims[[2]])
}


save_pheatmap_svg <- function(x, save_path, size = NA, width = NA, height = NA) {
    size <- decide_size(size = size[[1]], width = width, height = height)
    svg(save_path, width = size[[1]], height = size[[2]])
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}