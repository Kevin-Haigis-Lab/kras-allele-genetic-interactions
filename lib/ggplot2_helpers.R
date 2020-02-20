
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
    dims <- c()
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
    dims <- decide_size(size = size[[1]], width = width, height = height)

    if (tools::file_ext(save_path) != "svg") {
        message("File name does not have a 'svg' extension.")
    }

    ggsave(filename = save_path, plot = p,
           width = dims[[1]], height = dims[[2]])
}


#' A function factory for getting integer y-axis values.
integer_breaks <- function(n = 5, ..., rm_vals = c()) {
    fxn <- function(x) {
        breaks <- floor(pretty(x, n, ...))
        breaks <- breaks[!breaks %in% rm_vals]
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}


#' A function factory for getting rounded y-axis values.
round_breaks <- function(n = 5, n_dec = 2, ...) {
    fxn <- function(x) {
        breaks <- round(pretty(x, n, ...), n_dec)
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}


#' A function for `coord_trans()` to prevent `log(0)`.
#' use:
#'   ggplot(...) +
#'     ... +
#'     coord_trans(y = my_trans_log10) +
#'     .... +
my_trans_log10 <- scales::trans_new(
    name = "log10 pseudo-count +1",
    transform = function(x) { log10(x + 1) },
    inverse = function(x) { exp(x) - 1 }
)





#### ---- Saving a 'pheatmap' heatmap ---- ####

#' Save a pheatmap object
save_pheatmap_svg <- function(x, save_path, size = NA, width = NA, height = NA) {
    size <- decide_size(size = size[[1]], width = width, height = height)
    svg(save_path, width = size[[1]], height = size[[2]])
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

save_pheatmap_pdf <- function(x, save_path, size = NA, width = NA, height = NA) {
    size <- decide_size(size = size[[1]], width = width, height = height)
    pdf(save_path, width = size[[1]], height = size[[2]])
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}


#### ---- Adjust colors ---- ####

darken <- function(color, factor = 1.4){
    (col2rgb(color) / factor) %>%
        t() %>%
        rgb(maxColorValue = 255)
}

lighten <- function(colors, factor = 1.4){
    f <- function(color) {
        apply(col2rgb(color) * factor, 1, function(x) min(255, x)) %>%
        as.matrix() %>%
        t() %>%
        rgb(maxColorValue = 255)
    }
    purrr::map(colors, f)
}
