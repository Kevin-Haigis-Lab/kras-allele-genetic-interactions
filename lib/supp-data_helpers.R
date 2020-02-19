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
# section: section of the paper the table pertains to.
# number: order of the sheet in the section.
# sheet_name: the name of the final sheet in the Supp. Data file.
#
# Invisibly returns the original table.
save_supp_data <- function(tbl, section, number, sheet_name, verbose = TRUE) {

    file_name <- supp_data_filename(section, number, sheet_name)
    file_name <- sheets_path(file_name)

    janitor::clean_names(tbl) %>%
        write_tsv(file_name)

    if (verbose) {
        cat("Saved table for Supp. Data.\n")
        cat("sheet name: \"", clean_sheet_name(sheet_name), "\"\n", sep = "")
    }

    invisible(tbl)
}


# Get the file name for a sheet.
supp_data_filename <- function(section, number, sheet_name) {
    new_section <- str_pad(as.character(section), 2, "left", "0")
    new_number <- str_pad(as.character(number), 3, "left", "0")
    new_sheet_name <- clean_sheet_name(sheet_name)

    if (str_length(new_section) > 2) stop("Section value is too long.")
    if (str_length(new_number) > 3) stop("Number of sheet is too long.")
    if (str_length(new_sheet_name) < 3) stop("Sheet name is too short.")
    fname <- glue("{new_section}_{new_number}_{new_sheet_name}.txt")
    return(as.character(fname))
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


# Compile all of the separate sheets into a single Excel spreadsheet.
compile_supp_data <- function(verbose = TRUE) {
    all_sheet_paths <- list_all_sheet_paths()

    if (verbose) cat("Clearing old spreadsheet (if necessary).\n")
    clear_final_data_file()

    cat("Writing:\n")
    for (sheet_path in all_sheet_paths) {
        df <- read_tsv(sheet_path, col_types = cols())

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
    }
    if (verbose) cat("final file: \"", SUPP_DATA_FILE, "\"\n", sep = "")
}


# Return all sheet paths.
list_all_sheet_paths <- function() {
    list.files(SUPP_DATA_SHEETS_DIR, pattern = "txt$", full.names = TRUE)
}


# If the final spreadsheet already exists, remove it.
clear_final_data_file <- function() {
    if (file.exists(SUPP_DATA_FILE)) file.remove(SUPP_DATA_FILE)
}


# Prepare the sheet name from the file path.
prep_sheet_name <- function(sheet_path) {
    sheet_name <- file_sans_ext(basename(sheet_path))

    if (str_count(sheet_name, "_") != 2) stop("Sheet name does not conform.")

    sheet_name <- unlist(str_split_fixed(sheet_name, "_", 3)[, 3])
    str_replace_all(sheet_name, "-", " ")
}
