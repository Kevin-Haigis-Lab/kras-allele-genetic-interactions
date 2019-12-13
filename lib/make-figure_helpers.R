

build_comutation_figures <- function(numbers = NULL) {
    stop("This function has not been implemented.")
}

#### ---- Figure dimensions ---- ####

# Data sourced from NRJ guidelines.
FIG_WIDTHS <- c(88, 180)
FIG_HEIGHTS <- tibble::tribble(
    ~n_col, ~height_class, ~height, ~max_cation_wc,
         1,       "short",     130,            300,
         1,      "medium",     180,            150,
         1,        "tall",     220,             50,
         2,       "short",     185,            300,
         2,      "medium",     210,            150,
         2,        "tall",     225,             50
) %>%
    mutate(
        height_class_short = str_sub(height_class, 1, 1),
        height_class_num = rep(seq(1, 3), 2)
    )

#' Returns the figure dimensions (width and height) in "mm".
get_figure_dimensions <- function(n_col = 1, height_class = "medium") {

    width <- FIG_WIDTHS[[n_col]]
    height <- FIG_HEIGHTS %>%
        filter(n_col == !!n_col) %>%
        filter(
            height_class == !!height_class |
            height_class_short == !!height_class |
            height_class_num == !!height_class
        )

    if (nrow(height) == 0) {
        stop(glue("Height class '{height_class}' is not supported."))
    }

    return(list(
        width = width,
        height = height$height
    ))
}



#### ---- Saving and retrieving figure protos ---- ####

# The root directory for the protos.
FIGURE_DIR <- file.path("paper", "figures")
FIGURE_PROTOS_DIR <- file.path(FIGURE_DIR, "figure_protos")


#' Returns the numeric value for the lowest subversion of Figure `figure_num`.
get_latest_version_of_figure <- function(figure_num) {
    all_dirs <- list.dirs(FIGURE_PROTOS_DIR,
                          full.names = FALSE,
                          recursive = FALSE)
    idx <- str_detect(all_dirs, glue("_{figure_num}-"))
    all_dirs <- all_dirs[idx]
    subversions <- str_split_fixed(all_dirs, "_|-", 3)[, 3] %>%
        unlist() %>%
        as.numeric()
    return(max(subversions))
}


#' Get the path to proto file for version `version` of Figure `figure_num`.
#' The general format is "figure_protos/figure_00-000/filename.rds".
get_fig_proto_path <- function(name, figure_num, version = "latest") {

    base_n <- str_pad(as.character(figure_num), 2, pad = "0")
    if (version == "latest") {
        version <- get_latest_version_of_figure(base_n)
    }
    sub_n <- str_pad(as.character(version), 3, pad = "0")


    fig_dir <- glue("figure_{base_n}-{sub_n}")
    full_name <- glue("{name}.rds")

    file.path(FIGURE_PROTOS_DIR, fig_dir, full_name)
}


#' Get the path for the final figure file.
get_figure_path <- function(figure_num, version = "latest") {
    base_n <- str_pad(as.character(figure_num), 2, pad = "0")
    if (version == "latest") {
        version <- get_latest_version_of_figure(base_n)
    }
    sub_n <- str_pad(as.character(version), 3, pad = "0")

    fig_dir <- glue("figure_{base_n}-{sub_n}")
    versioned_name <- glue("Figure_{base_n}-{sub_n}.svg")
    unversioned_name <- glue("Figure_{base_n}.svg")

    return(list(
        versioned = file.path(FIGURE_PROTOS_DIR, fig_dir, versioned_name),
        unversioned = file.path(FIGURE_DIR, unversioned_name)
    ))
}


#' Save the final SVG for Figure `figure_num`.
#' Two files are saved, one with the version number and one without.
save_figure <- function(p,
                        figure_num, version = "latest",
                        dim=NULL,
                        n_col=NULL, height_class = NULL,
                        unversioned_only = FALSE) {
    if (is.null(dim)) {
        dim <- get_figure_dimensions(n_col, height_class)
    }

    file_names <- get_figure_path(figure_num, version)

    # Save versioned.
    ggsave(
        file_names$versioned, p,
        width = dim$width, height = dim$height, unit = "mm"
    )
    if (!unversioned_only) {
        ggsave(
            file_names$unversioned, p,
            width = dim$width, height = dim$height, unit = "mm"
        )
    }
}


#' Read in the proto RDS file names `name` for version `version` of
#' Figure `figure_num`.
read_fig_proto <- function(name, figure_num, version = "latest") {
    path <- get_fig_proto_path(name, figure_num, version)
    readRDS(path)
}
