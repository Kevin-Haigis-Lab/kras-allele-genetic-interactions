# Preparation of input data for GSEA

# Paths to save the expression and cls files to.
gsea_save_paths <- function(cancer, allele) {
    parent_dir <- file.path("data", "gsea", "input")
    file_name <- glue("{cancer}_{allele}")
    return(list(
        expr = file.path(parent_dir, paste0(file_name, ".txt")),
        cls = file.path(parent_dir, paste0(file_name, ".cls"))
    ))
}



# Prepare the expression input data for GSEA in the following format.
#    Missing values are turned into blanks ("").
# | NAME  | DESCRIPTION | sample1 | sample2 | sample3 |
# |-------|-------------|---------|---------|---------|
# | gene1 | na          |         |         |         |
# | gene2 | na          |         |         |         |
make_expression_dataset <- function(data) {
    wide_data <- data %>%
        select(dep_map_id, hugo_symbol, gene_effect) %>%
        unique() %>%
        group_by(dep_map_id, hugo_symbol) %>%
        filter(n() == 1) %>%
        ungroup() %>%
        pivot_wider(
            names_from = dep_map_id,
            values_from = gene_effect
        ) %>%
        add_column(DESCRIPTION = "na") %>%
        dplyr::rename(NAME = "hugo_symbol") %>%
        select(NAME, DESCRIPTION, everything())
    wide_data[is.na(wide_data)] <- ""
    invisible(wide_data)
}

# Make the CLS file for GSEA.
make_cls_file <- function(wide_data, data, output_path) {
    # collect line 1 info
    number_of_samples <- ncol(wide_data) - 2
    number_of_classes <- 2

    # collect line 2 info
    classes <- unique(data$allele)
    class_0_name <- classes[classes != "other"]
    class_1_name <- "other"

    # collect line 3 info
    cls_nums <- tibble(dep_map_id = colnames(wide_data)[c(-1, -2)]) %>%
        left_join(
            { data %>% select(dep_map_id, allele) %>% unique() },
            by = "dep_map_id"
        ) %>%
        mutate(
            cls_num = ifelse(allele == class_0_name, 0, 1)
        ) %>%
        pull(cls_num)

    # write out
    if (file.exists(output_path)) { file.remove(output_path) }
    cat(number_of_samples, number_of_classes, "1\n", file = output_path)
    cat("#", class_0_name, class_1_name, "\n", file = output_path, append = TRUE)
    cat(cls_nums, "\n", file = output_path, append = TRUE)
}

# Make the input files for GSEA for a cancer and allele
make_gsea_input_files <- function(cancer, allele, data) {
    # destinations for expression and CLS files
    save_paths <- gsea_save_paths(cancer, allele)

    # data for the cancer and adjusted allele names to <allele> and "other"
    mod_data <- data %>%
        filter(cancer == !!cancer) %>%
        mutate(allele = ifelse(allele == !!allele, allele, "other"))

    # expression file
    expr_df <- make_expression_dataset(mod_data)
    write_tsv(expr_df, save_paths$expr)

    # CLS file
    make_cls_file(expr_df, mod_data, save_paths$cls)

    return(NULL)
}

allele_cancer_tib <- model_data %>%
    filter(cancer != "MM") %>%
    select(cancer, allele) %>%
    unique() %>%
    pwalk(make_gsea_input_files, data = model_data)



# model_data %>%
#     select(hugo_symbol) %>%
#     mutate(
#         `Probe Set ID` = hugo_symbol,
#         `Gene Symbol` = hugo_symbol,
#         `Gene Title` = "some text"
#     ) %>%
#     select(-hugo_symbol) %>%
#     write_tsv(file.path("data", "gsea", "hugo_symbol.chip"))
