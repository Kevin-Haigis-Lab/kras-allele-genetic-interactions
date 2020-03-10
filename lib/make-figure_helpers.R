
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


#### ---- Figure base theme ---- ####

theme_comutation <- function() {
    theme_bw(base_size = 6, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5, hjust = 1),
        axis.text.x = element_text(size = 5, vjust = 1),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-2, -2, -2, -2, "mm")),
        strip.background = element_blank()
    )
}


theme_classic_comutation <- function() {
    theme_classic(base_size = 6, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5, hjust = 1),
        axis.text.x = element_text(size = 5, vjust = 1),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-2, -2, -2, -2, "mm")),
        strip.background = element_blank()
    )
}

theme_minimal_comutation <- function() {
    theme_minimal(base_size = 6, base_family = "Arial") %+replace%
    theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5, hjust = 1),
        axis.text.x = element_text(size = 5, vjust = 1),
        axis.ticks = element_blank(),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-2, -2, -2, -2, "mm")),
        strip.background = element_blank()
    )
}

theme_graph_comutation <- function() {
    theme_graph(base_size = 6, base_family = "Arial") %+replace%
    theme(
        plot.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        plot.tag = element_text(size = 7,
                                face = "bold",
                                margin = margin(-2, -2, -2, -2, "mm"))
    )
}



#### ---- Saving and retrieving figure protos ---- ####

# The root directory for the protos.
FIGURE_DIR <- file.path("paper", "figures")
FIGURE_PROTOS_DIR <- file.path(FIGURE_DIR, "figure_protos")


#' Returns the numeric value for the lowest subversion of Figure `figure_num`.
get_latest_version_of_figure <- function(figure_num, supp) {
    all_dirs <- list.dirs(FIGURE_PROTOS_DIR,
                          full.names = FALSE,
                          recursive = FALSE)
    fig_type <- ifelse(supp, "suppfigure", "figure")
    idx <- str_detect(all_dirs, glue("{fig_type}_{figure_num}-"))
    all_dirs <- all_dirs[idx]
    subversions <- str_split_fixed(all_dirs, "_|-", 3)[, 3] %>%
        unlist() %>%
        as.numeric()
    return(max(subversions))
}


#' Get the path to proto file for version `version` of Figure `figure_num`.
#' The general format is "figure_protos/figure_00-000/filename.rds".
get_fig_proto_path <- function(name,
                               figure_num,
                               version = "latest",
                               supp = FALSE) {
    name <- basename(name)
    name <- file_sans_ext(name)

    base_n <- str_pad(as.character(figure_num), 2, pad = "0")
    if (version == "latest") {
        version <- get_latest_version_of_figure(base_n, supp = supp)
    }
    sub_n <- str_pad(as.character(version), 3, pad = "0")

    if (!supp) {
        fig_dir <- glue("figure_{base_n}-{sub_n}")
        full_name <- glue("{name}.rds")
    } else {
        fig_dir <- glue("suppfigure_{base_n}-{sub_n}")
        full_name <- glue("{name}.rds")
    }


    file.path(FIGURE_PROTOS_DIR, fig_dir, full_name)
}


#' Get the path for the final figure file.
get_figure_path <- function(figure_num,
                            version = "latest",
                            supp = FALSE,
                            file_fmt = "svg",
                            subdir = "") {
    base_n <- str_pad(as.character(figure_num), 2, pad = "0")
    if (version == "latest") {
        version <- get_latest_version_of_figure(base_n, supp = supp)
    }
    sub_n <- str_pad(as.character(version), 3, pad = "0")

    if (!supp) {
        fig_dir <- glue("figure_{base_n}-{sub_n}")
        versioned_name <- glue("Figure_{base_n}-{sub_n}.{file_fmt}")
        unversioned_name <- glue("Figure_{base_n}.{file_fmt}")
    } else {
        fig_dir <- glue("suppfigure_{base_n}-{sub_n}")
        versioned_name <- glue("SuppFigure_{base_n}-{sub_n}.{file_fmt}")
        unversioned_name <- glue("SuppFigure_{base_n}.{file_fmt}")
    }

    return(list(
        versioned = file.path(FIGURE_PROTOS_DIR, fig_dir, versioned_name),
        unversioned = file.path(FIGURE_DIR, subdir, unversioned_name)
    ))
}


#' Save the final SVG for Figure `figure_num`.
#' Two files are saved, one with the version number and one without.
save_figure <- function(p,
                        figure_num,
                        version = "latest",
                        supp = FALSE,
                        dim = NULL,
                        n_col = NULL,
                        height_class = NULL,
                        unversioned_only = FALSE) {
    if (is.null(dim)) {
        dim <- get_figure_dimensions(n_col, height_class)
    }

    svg_names <- get_figure_path(figure_num,
                                 version,
                                 supp = supp,
                                 file_fmt = "svg")

    jpg_names <- get_figure_path(figure_num,
                                 version,
                                 supp = supp,
                                 file_fmt = "jpeg")

    pdf_names <- get_figure_path(figure_num,
                                 version,
                                 supp = supp,
                                 file_fmt = "pdf",
                                 subdir = "pdfs")

    # Save versioned.
    for (names_list in list(svg_names, jpg_names)) {
        ggsave(
            names_list$versioned, p,
            width = dim$width, height = dim$height, unit = "mm"
        )
    }

    if (!unversioned_only) {
        for (names_list in list(svg_names, jpg_names, pdf_names)) {
            ggsave(
                names_list$unversioned, p,
                width = dim$width, height = dim$height, unit = "mm"
            )
        }
    }

    for (name in pdf_names$unversioned) {
        ggsave(
            name, p,
            width = dim$width, height = dim$height, unit = "mm",
            device = cairo_pdf
        )
    }

    # Save a JPEG to "reports/content/home/gallery/gallery/"
    gallery_path <- file.path(
        "reports", "content", "home", "gallery", "gallery",
        basename(jpg_names$unversioned)
    )
    ggsave(
        gallery_path, p,
        width = dim$width, height = dim$height, unit = "mm"
    )

}


#' Read in the proto RDS file names `name` for version `version` of
#' Figure `figure_num`.
read_fig_proto <- function(name, figure_num, version = "latest", supp = FALSE) {
    path <- get_fig_proto_path(name, figure_num, version, supp = supp)
    readRDS(path)
}


#### ---- Start a new figure ---- ####


TEMPLATE_PATH <- file.path(FIGURE_PROTOS_DIR, "_make-figure-template.R")

glue_figure_template <- function(figure_num, version, supp) {
    # Make `figure name` and `theme_suffix`
    if (supp) {
        figure_name <- glue("Supplemental Figure {figure_num}")
        theme_suffix <- glue("S{figure_num}")
    } else {
        figure_name <- glue("Figure {figure_num}")
        theme_suffix <- as.character(figure_num)
    }
    version <- as.character(version)
    supp_chr <- ifelse(supp, "TRUE", "FALSE")
    figure_num <- as.character(figure_num)

    template <- readLines(TEMPLATE_PATH)
    for (i in seq(1, length(template))) {

        cond1 <- str_length(template[[i]]) == 0
        cond2 <- !str_detect(template[[i]], "__")
        if (cond1 | cond2) next

        template[[i]] <- template[[i]] %>%
            str_replace("__FIGURE_NAME__", figure_name) %>%
            str_replace("__FIGURE_NUM__", figure_num) %>%
            str_replace("__SUPP__", supp_chr) %>%
            str_replace("__VERSION__", as.character(version)) %>%
            str_replace("__THEME_SUFFIX__", theme_suffix)
    }
    return(template)
}


initialize_figure <- function(figure_num,
                              version = "next",
                              supp = FALSE) {

    # Get version number
    if (version == "next") {
        version <- suppressWarnings(get_latest_version_of_figure(figure_num,
                                                                 supp))
        version <- version + 1
        if (is.infinite(version)) version <- 1
    }

    # Get the correct strings for figure number and version.
    base_n <- str_pad(as.character(figure_num), 2, pad = "0")
    sub_n <- str_pad(as.character(version), 3, pad = "0")

    # Get path name depending on if it is supplementary or not.
    if (supp) {
        file_name <- glue("make-suppfigure_{base_n}-{sub_n}.R")
    } else {
        file_name <- glue("make-figure_{base_n}-{sub_n}.R")
    }
    path <- file.path(FIGURE_PROTOS_DIR, file_name)

    # Stop if the file already exists.
    # Do not overwrite already existing files.
    if (file.exists(path)) {
        stop(glue("File already exists for Fig. {figure_num}, v. {sub_n}"))
    }

    # Make directory for the protos for the figure.
    fig_dir <- str_remove(path, "make-")
    fig_dir <- str_remove(fig_dir, "\\.R$")
    if (!dir.exists(fig_dir)) dir.create(fig_dir)

    # Insert the correct information into the template.
    template <- glue_figure_template(figure_num, version, supp)
    writeLines(template, path)
}


#### ---- Build all figures ---- ####

#' Extract the title, number, and version of the make script.
#' Information is returned as a labeled list.
get_figure_info_from_name <- function(make_file_name) {
    info <- list()
    fn <- basename(file_sans_ext(make_file_name))
    info$title <- ifelse(str_detect(fn, "supp"), "Supp. Figure", "Figure")
    info$number <- as.numeric(str_extract(fn, "(?<=figure_)[:digit:]+(?=-)"))
    info$version <- as.numeric(str_extract(fn, "(?<=\\-)[:digit:]+$"))
    return(info)
}


# Build all of the figures.
build_comutation_figures <- function() {
    figure_make_scripts <- list.files(FIGURE_PROTOS_DIR, full.names = TRUE)
    idx <- str_detect(basename(figure_make_scripts), "^make.*\\.R$")
    figure_make_scripts <- figure_make_scripts[idx]

    for (f in figure_make_scripts) {
        info <- get_figure_info_from_name(f)
        msg <- glue("=> Making {info$title} {info$number} (v{info$version})")
        cat(str_pad("", str_length(msg), pad = "="), "\n")
        cat(msg, "\n")

        source(f)

        cat(str_pad("", str_length(msg), pad = "="), "\n\n")
    }
    cat("Done\n")
}
