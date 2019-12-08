
# The root directory for the protos.
FIGURE_PROTOS_DIR <- file.path("paper", "figures", "figure_protos")

#' Get the path to proto file for version `version` of Figure `figure_num`.
get_fig_proto_path <- function(name, figure_num, version = 1) {
    fig_dir <- glue("figure_{figure_num}-{figure_num}")
    full_name <- glue("{name}.rds")
    file.path(FIGURE_PROTOS_DIR, fig_dir, name)
}

