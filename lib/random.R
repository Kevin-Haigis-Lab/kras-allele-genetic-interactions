
# Count the number of lines with text in all the files in a directory.
line_count_per_directory <- function(directory) {
    counter <- 0
    for (file in list.files(directory, full.names = TRUE, pattern = "R")) {
        all_lines <- readLines(file, warn = FALSE) %>% unlist() %>% str_trim()
        all_lines <- all_lines[str_length(all_lines) > 0]
        counter <- counter + length(all_lines)
    }
    return(counter)
}



# The number of lines of code in this project
comutation_proj_lines_of_code <- function() {
    lapply(
            c("src", "munge", "lib", "tests"),
            line_count_per_directory
        ) %>%
        unlist() %>%
        sum()
}