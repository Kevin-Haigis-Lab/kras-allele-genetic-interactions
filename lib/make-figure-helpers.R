
#### ---- Figure dimensions ---- ####

# Data sourced from NRJ guidelines.
FIG_WIDTHS <- c(88, 180)
FIG_HEIGHTS <- tibble::tribble(
  ~n_col, ~height_class, ~height, ~max_cation_wc,
  1, "short", 130, 300,
  1, "medium", 180, 150,
  1, "tall", 220, 50,
  2, "short", 185, 300,
  2, "medium", 210, 150,
  2, "tall", 225, 50
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
      axis.title = element_text(size = 7),
      axis.text.y = element_text(size = 6, hjust = 1),
      axis.text.x = element_text(size = 6, vjust = 1),
      axis.ticks = element_blank(),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(-2, -2, -2, -2, "mm")
      ),
      strip.background = element_blank()
    )
}


theme_classic_comutation <- function() {
  theme_classic(base_size = 6, base_family = "Arial") %+replace%
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text.y = element_text(size = 6, hjust = 1),
      axis.text.x = element_text(size = 6, vjust = 1),
      axis.ticks = element_blank(),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(-2, -2, -2, -2, "mm")
      ),
      strip.background = element_blank()
    )
}

theme_minimal_comutation <- function() {
  theme_minimal(base_size = 6, base_family = "Arial") %+replace%
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text.y = element_text(size = 6, hjust = 1),
      axis.text.x = element_text(size = 6, vjust = 1),
      axis.ticks = element_blank(),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(-2, -2, -2, -2, "mm")
      ),
      strip.background = element_blank()
    )
}

theme_graph_comutation <- function() {
  theme_graph(base_size = 6, base_family = "Arial") %+replace%
    theme(
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      plot.tag = element_text(
        size = 7,
        face = "bold",
        margin = margin(-2, -2, -2, -2, "mm")
      )
    )
}



#### ---- Saving and retrieving figure protos ---- ####

# Key directories
FIGURE_DIR <- file.path("paper", "figures")
FIGURE_PROTOS_DIR <- file.path(FIGURE_DIR, "make_figure_scripts")
FIGURE_IMAGE_DIR <- file.path(FIGURE_DIR, "all_images")
GGPROTO_OBJECTS_DIR <- file.path(FIGURE_DIR, "ggproto_objects")


# Save a ggproto object to the directory for the figures.
saveFigRds <- function(obj, name) {
  name <- file_sans_ext(basename(name))
  saveRDS(obj, file.path(GGPROTO_OBJECTS_DIR, name))
}

# Get he figure number padded to size 3.
pad_fig_number <- function(num) {
  str_pad(as.character(num), 3, side = "left", pad = "0")
}


# Get the path to the make figure file for a number.
get_figure_file_path <- function(num) {
  if (num > 999) {
    stop("No file can be numbered above 999.")
  }
  num_pad <- pad_fig_number(num)
  file.path(FIGURE_PROTOS_DIR, glue("make-figure_{num_pad}.R"))
}


# Returns a vector of all of the figure files.
get_all_figure_files <- function(full_names = TRUE) {
  current_files <- list.files(FIGURE_PROTOS_DIR, full.names = full_names)
  current_files <- unlist(current_files)
  idx <- str_detect(current_files, "make-figure_[:digit:]{3}\\.R")
  current_files <- current_files[idx]
}


# Get the next figure number.
get_max_figure_number <- function() {
  current_files <- get_all_figure_files(full_names = FALSE)

  # First file.
  if (length(current_files) == 0) {
    return(0)
  }

  current_num <- str_extract(current_files, "[:digit:]{3}")
  current_num <- as.numeric(current_num)
  return(max(current_num))
}


# Get the output path for a figure image.
get_figure_img_path <- function(num, file_fmt) {
  file.path(
    FIGURE_IMAGE_DIR,
    glue("figure_{pad_fig_number(num)}.{file_fmt}")
  )
}

# Save the final SVG for Figure `figure_num`.
# Two files are saved, one with the version number and one without.
save_figure <- function(p,
                        figure_num,
                        dim = NULL,
                        n_col = NULL,
                        height_class = NULL) {
  if (is.null(dim)) {
    dim <- get_figure_dimensions(n_col, height_class)
  }

  svg_name <- get_figure_img_path(figure_num, file_fmt = "svg")
  jpg_name <- get_figure_img_path(figure_num, file_fmt = "jpeg")

  for (f in c(svg_name, jpg_name)) {
    ggsave(
      filename = f,
      plot = p,
      width = dim$width, height = dim$height,
      unit = "mm", dpi = 500
    )
  }
  invisible(NULL)
}


# Read in the proto RDS file names `name`.
read_fig_proto <- function(name) {
  name <- file_sans_ext(basename(name))
  path <- file.path(GGPROTO_OBJECTS_DIR, name)
  readRDS(path)
}


#### ---- Start a new figure ---- ####


TEMPLATE_PATH <- file.path(FIGURE_PROTOS_DIR, "_make-figure-template.R.template")

glue_figure_template <- function(figure_num) {
  figure_num_pad <- pad_fig_number(figure_num)
  figure_name <- glue("Figure {figure_num_pad}")
  figure_num <- as.character(figure_num)

  template <- readLines(TEMPLATE_PATH)
  for (i in seq(1, length(template))) {
    cond1 <- str_length(template[[i]]) == 0
    cond2 <- !str_detect(template[[i]], "__")
    if (cond1 | cond2) next

    template[[i]] <- template[[i]] %>%
      str_replace("__FIGURE_NAME__", figure_name) %>%
      str_replace("__FIGURE_NUM__", figure_num) %>%
      str_replace("__THEME_SUFFIX__", figure_num)
  }
  return(template)
}


initialize_figure <- function() {
  next_num <- get_max_figure_number() + 1
  file_name <- get_figure_file_path(next_num)

  # Stop if the file already exists.
  # Do not overwrite already existing files.
  if (file.exists(file_name)) {
    stop(glue(
      "Something has gone wrong - file for {next_num} already exists!"
    ))
  } else {
    message(glue("Creating {basename(file_name)}"))
  }

  # Insert the correct information into the template.
  template <- glue_figure_template(next_num)
  writeLines(template, file_name)
}


#### ---- Build all figures ---- ####


source_comutation_figure <- function(f) {
  fig_num <- str_extract(basename(f), "[:digit:]{3}")
  msg <- glue("=> Making Figure {fig_num}")
  cat(str_rep("=", str_length(msg)), "\n")
  cat(msg, "\n")

  source(f)

  cat(str_rep("=", str_length(msg)), "\n\n")
}

# Build one or more specific figures.
build_comutation_figure <- function(nums) {
  if (is.infinite(nums)) {
    build_comutation_figures()
  } else {
    for (i in nums) {
      f <- get_figure_file_path(i)
      if (file.exists(f)) {
        source_comutation_figure(f)
      } else {
        cat(glue("(Skipping {basename(f)} - does not exist)"), "\n")
      }
    }
  }
}


# Build all of the figures.
build_comutation_figures <- function() {
  clear_source_data()
  figure_make_scripts <- get_all_figure_files()

  for (f in figure_make_scripts) {
    source_comutation_figure(f)
  }

  cat("Done\n")
}


#### ---- Compiling final figures ---- ####

JSON_PATH <- file.path(FIGURE_DIR, "figure_conversion.json")

check_conversion_df <- function(df) {
  stopifnot(nrow(df) == nrow(unique(df)))
  stopifnot(n_distinct(df$make_num) == nrow(df))
  return(df)
}


get_conversion_df <- function() {
  jsonlite::fromJSON(JSON_PATH) %>%
    as_tibble() %>%
    check_conversion_df()
}


get_final_figure_path <- function(num, supp, fmt = "svg") {
  num <- str_pad(num, width = 2, side = "left", pad = "0")
  prefix <- ifelse(supp, "Supp_Fig", "Fig")
  file.path(FIGURE_DIR, glue("{prefix}_{as.character(num)}.{fmt}"))
}


copy_final_figure <- function(figure_num, make_num, supp) {
  for (fmt in c("svg", "jpeg")) {
    base_img_path <- get_figure_img_path(make_num, fmt)
    output_path <- get_final_figure_path(figure_num, supp, fmt)
    file.copy(
      base_img_path,
      output_path,
      overwrite = TRUE,
      recursive = FALSE
    )
  }
  invisible(NULL)
}


make_final_pdfs <- function(...) {
  # system("paper/figures_to_pdf.sh")
  message("Note that PDFs must be made locally due to font issues.")
}


remove_old_final_figures <- function() {
  old_figure_files <- unlist(list.files(FIGURE_DIR, full.names = TRUE))

  idx <- str_detect(old_figure_files, "Fig_[:digit:]+\\.(jpeg|pdf|svg)")
  old_figure_files <- old_figure_files[idx]

  a <- file.remove(old_figure_files)
  invisible(NULL)
}

copy_final_figures <- function() {
  remove_old_final_figures()
  pwalk(get_conversion_df(), copy_final_figure)
  make_final_pdfs()
}
