# The API for saving data to Supplementary Data and compiling the final
# Excel spreadsheet.


#### ---- Folder structure ---- ####

SUPP_DATA_DIR <- file.path("paper", "supplemental_data")
SUPP_DATA_SHEETS_DIR <- file.path(SUPP_DATA_DIR, "sheets")

for (a in c(SUPP_DATA_DIR, SUPP_DATA_SHEETS_DIR)) {
    if (!dir.exists(a)) dir.create(a)
}


# Append the path the the "supplemental_data" directory to `...`.
suppdata_path <- function(...) {
    file.path(SUPP_DATA_DIR, ...)
}


# Append the path the the "sheets" subdirectory to `...`.
sheets_path <- function(...) {
    file.path(SUPP_DATA_SHEETS_DIR, ...)
}


#### ---- Writing data frames to sheets ---- ####

# Write a table to be included as a sheet in the Supp. Data.
#
# tbl: a data frame to be written as a TSV.
# num: a number to assign to the sheet - ordering can be changed later.
# sheet_name: the name of the final sheet in the Supp. Data file.
#
# Invisibly returns the original table.
save_supp_data <- function(tbl, num, sheet_name, verbose = TRUE) {

    file_name <- supp_data_filename(num, sheet_name)
    file_name <- sheets_path(file_name)

    janitor::clean_names(tbl) %>%
        write_tsv(file_name)

    if (verbose) {
        cat("Saved table for Supp. Data.\n")
        cat("sheet name: \"", clean_sheet_name(sheet_name), "\"\n", sep = "")
    }

    invisible(tbl)
}


# Returns the maximum number of the currently saved sheets.
get_max_sheet_number <- function() {
    all_files <- list_all_sheet_paths()
    if (length(all_files) == 0) {
        return(0)
    } else {
        all_nums <- extract_sheet_nums(all_files) %>%
            as.numeric()
        return(max(all_nums))
    }
}


# Get the file name for a sheet.
supp_data_filename <- function(number, sheet_name) {
    new_number <- str_pad(as.character(number), 3, "left", "0")
    new_sheet_name <- clean_sheet_name(sheet_name)

    if (str_length(new_number) > 3) stop("Number of sheet is too long.")
    if (str_length(new_sheet_name) < 3) stop("Sheet name is too short.")
    fname <- glue("{new_number}_{new_sheet_name}.txt")
    return(as.character(fname))
}


# Get the sheet numbers of the files.
extract_sheet_nums <- function(fs) {
    unlist(fs) %>%
        basename() %>%
        str_extract("^[:digit:]{3}")
}


# Clean the input sheet name for the file name.
clean_sheet_name <- function(sheet_name) {
    clean_name <- str_to_lower(sheet_name)
    for (s in c(" ", "_")) {
        clean_name <- str_replace_all(clean_name, s, "-")
    }
    return(clean_name)
}


#### ---- Compiling sheet ---- ####

# The final Excel spreadsheet.
SUPP_DATA_FILE <- file.path(SUPP_DATA_DIR, "Supplemental-Data.xlsx")
SUPP_DATA_ORDER_JSON <- file.path(SUPP_DATA_DIR, "suppdata-sheet-order.json")


# Compile all of the separate sheets into a single Excel spreadsheet.
compile_supp_data <- function(verbose = TRUE) {
    all_sheet_paths <- list_all_sheet_paths() %>% put_sheets_in_order()

    if (verbose) cat("Clearing old spreadsheet (if necessary).\n")
    clear_final_data_file()

    if (verbose) cat("Apportioning memory for java.\n")
    options(java.parameters = "-Xmx8g")
    clear_java_memory()

    cat("Writing:\n")
    for (sheet_path in all_sheet_paths) {
        df <- read_tsv(sheet_path, col_types = cols(.default = "c"))

        sheet_name <- prep_sheet_name(sheet_path)
        if (verbose) cat("  -->", sheet_name, "\n")

        xlsx::write.xlsx(
            x = df,
            file = SUPP_DATA_FILE,
            sheetName = prep_sheet_name(sheet_path),
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE
        )
        clear_java_memory()
    }
    if (verbose) cat("final file: \"", SUPP_DATA_FILE, "\"\n", sep = "")
}


clear_java_memory <- function() {
    gc()
    XLConnect::xlcFreeMemory()
    invisible(NULL)
}


# Get a data frame mapping the sheet number to the final number.
get_sheet_order_df <- function() {
    jsonlite::fromJSON(SUPP_DATA_ORDER_JSON) %>%
        as_tibble()
}


# Get a tibble of the supp. data sheet names and their make numbers.
get_sheet_file_numbers <- function(fs) {
    tibble(sheet = fs) %>%
        mutate(make_num = extract_sheet_nums(fs),
               make_num = as.numeric(make_num))
}


# Put the sheets in the correct order using `SUPP_DATA_ORDER_JSON`.
put_sheets_in_order <- function(fs) {
    left_join(get_sheet_file_numbers(fs),
              get_sheet_order_df(),
              by = "make_num") %>%
        arrange(final_order) %>%
        pull(sheet)
}


# Return all sheet paths in correct order.
# The order can be adjusted at `SUPP_DATA_ORDER_JSON`.
list_all_sheet_paths <- function() {
    list.files(SUPP_DATA_SHEETS_DIR,
               pattern = "txt$",
               full.names = TRUE)

}


# If the final spreadsheet already exists, remove it.
clear_final_data_file <- function() {
    if (file.exists(SUPP_DATA_FILE)) { file.remove(SUPP_DATA_FILE) }
    invisible(NULL)
}



replace_specific_terms <- function(s) {
    df <- tibble::tribble(
        ~ old, ~ new,
        " coad", " COAD",
        " luad", " LUAD",
        " mm", " MM",
        " paad", " PAAD",
        "kras", "KRAS"
    )

    for (i in seq(1, nrow(df))) {
        s <- str_replace_all(s, df$old[[i]], df$new[[i]])
    }

    return(s)
}


# Prepare the sheet name from the file path.
prep_sheet_name <- function(sheet_path) {
    sheet_name <- file_sans_ext(basename(sheet_path))
    sheet_name <- unlist(str_split_fixed(sheet_name, "_", 2)[, 2])
    str_replace_all(sheet_name, "-", " ") %>%
        replace_specific_terms()
}
