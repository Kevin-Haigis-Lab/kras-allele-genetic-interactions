
temp_sim_dir <- file.path("temp", "comut-fdr")
temp_input_dir <- file.path(temp_sim_dir, "input")
temp_output_dir <- file.path(temp_sim_dir, "output")

simulation_input_file_name <- function(idx) {
    file.path(temp_input_dir, glue("simulation_{idx}.tsv"))
}

simulation_output_file_name <- function(idx) {
    file.path(temp_output_dir, glue("simulation_{idx}.qs"))
}
