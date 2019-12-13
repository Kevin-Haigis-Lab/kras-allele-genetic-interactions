
#### ---- Globally relevant functions ---- ####

#' Get Hugo identifiers from the DepMap combined name format
#'
#' Extracts the Hugo gene identifiers from a vector of values with the format:
#' "Hugo (Entrez)".
#'
get_hugo_from_depmap_ids <- function(x, split_pattern = " ") {
    y <- unlist(stringr::str_split_fixed(x, split_pattern, n = 2)[, 1])
    return(y)
}


#' Get Entrez identifiers from the DepMap combined name format
#'
#' Extracts the Hugo gene identifiers from a vector of values with the format:
#' "Hugo (Entrez)".
#'
get_entrez_from_depmap_ids <- function(x, convert_to_num = FALSE) {
    y <- stringr::str_extract(x, "(?<=\\()[:digit:]+(?=\\))")
    if (convert_to_num) { y <- as.numeric(y) }
    return(y)
}


#' Standard error of the mean
#'
#' Calculate the standard error of the mean.
#'
sem <- function(x) sd(x) / sqrt(length(x))


#' Create or reset a directory in 'graphs'
#'
#' This is useful to call at the top of an analysis to remove the old graphs
#' such that there are no phantom results from old analyses.
reset_graph_directory <- function(dir_name) {
    dir_path <- file.path("graphs", dir_name)
    reset_directory(dir_path, "graph")
}


#' Create or reset a directory in 'tables'
#'
#' This is useful to call at the top of an analysis to remove the old tables
#' such that there are no phantom results from old analyses.
reset_table_directory <- function(dir_name) {
    dir_path <- file.path("tables", dir_name)
    reset_directory(dir_path, "table")
}


#' Create or reset a directory
#'
#' Do not call this directly. There are helpers for "graphs" and "tables".
#' If additional result locations need to be accessed, create a helper function
#' to maintain the API.
reset_directory <- function(dir_path, location) {
    message(glue("Reseting {location} directory: {basename(dir_path)}"))
    if (dir.exists(dir_path)) {
        unlink(dir_path, recursive = TRUE)
    }
    dir.create(dir_path)
}


#' Source the files in the 'lib' directory.
source_lib <- function() {
    for (f in list.files("lib", full.name = TRUE)) source(f)
}


#' Get the full path for a plot name.
plot_path <- function(...) {
    file.path("graphs", ...)
}


#' Get the full path for a plot name.
table_path <- function(...) {
    file.path("tables", ...)
}



#### ---- Helpful regular expressions ---- ####

# helpful data for identifying KRAS mutations
kras_hotspot_codons <- list(num = c(12, 13, 61, 117, 146))
kras_hotspot_codons$char <- as.character(kras_hotspot_codons$num)
kras_hotspot_codons$regex <- paste0(kras_hotspot_codons$char, collapse = "|")

# cancer regex
cancer_regex <- "COAD|LUAD|PAAD|MM|SKCM"



#### ---- Helpful colors/palettes ---- ####


# coding mutation variant classification terms
coding_mut_var_classes <- c(
    "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
    "Nonstop_Mutation", "Splice_Site", "Start_Codon_Del", "Start_Codon_Ins",
    "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins"
)

mapping_mutation_types_to_human_readable <- list(
    "missense_mutation" = "missense mutation",
    "frame_shift_ins" = "frameshift insertion",
    "frame_shift_del" = "frameshift deletion",
    "3'utr" = "3' UTR",
    "intron" = "intronic",
    "in_frame_del" = "in-frame deletion",
    "5'utr" = "5' UTR",
    "5'flank" = "5' flank",
    "silent" = "silent",
    "3'flank" = "3' flank",
    "splice_site" = "splice-site",
    "nonstop_mutation" = "nonstop mutation",
    "nonsense_mutation" = "nonsense mutation",
    "in_frame_ins" = "in-frame insertion",
    "translation_start_site" = "translation start site",
    "splice_region" = "splice-site",
    "rna" = "RNA",
    "igr" = "IGR",
    "intronic" = "intronic",
    "3flank" = "3' flank",
    "5flank" = "5' flank",
    "3utr" = "3' UTR",
    "flank" = "flank",
    "5utr" = "5' UTR",
    "inframe_insertion" = "in-frame insertion",
    "frameshift_insertion" = "frameshift insertion",
    "inframe_deletion" = "in-frame deletion",
    "frameshift_deletion" = "frameshift deletion",
    "unknown" = "unknown",
    "utr" = "UTR",
    "targeted_region" = "targeted region",
    "de_novo_start_outofframe" = "de novo start (out-of-frame)",
    "start_codon_snp" = "start codon SNV",
    "de_novo_start_inframe" = "de novo start (in-frame)",
    "missense" = "missense mutation"
)


## make color palette: does NOT run on O2, but does on Mac
# s <- 0
# while(TRUE) {
#     cat("Using", s, "as the seed value. ")
#     set.seed(s)
#     unique_muts <- unique(unname(unlist(mapping_mutation_types_to_human_readable)))
#     n <- length(unique_muts)
#     cols <- randomcoloR::randomColor(n, luminosity = "light")
#     scales::show_col(cols)
#     a <- readline(prompt = "Press Enter to continue... ")
#     if (a == "s") { break }
#     s <- s + 1
# }
#
# set.seed(30)
# unique_muts <- unique(unname(unlist(mapping_mutation_types_to_human_readable)))
# n <- length(unique_muts)
# cols <- randomcoloR::randomColor(n)
# for (i in 1:n) {
#     if (i == 1) { cat("\t") }

#     if (i %% 6 == 0) { cat("\n\t") }
#     cat("\"", cols[[i]], "\", ", sep = "")

#     if (i == n) { cat("\n") }
# }

mutation_pal <- list(
    "#fce094", "#bf7909", "#a4f9c8", "#8cf2b3", "#7e5cb5",
    "#dd8c2e", "#be2ecc", "#c254d8", "#e27bfc", "#1372d8", "#110bc4",
    "#eef9a7", "#075f82", "#fc6ab3", "#917add", "#fca69f", "#d67b0c",
    "#a8eae9", "#d8841e", "#e51c02", "#d4b9f7", "#d87b34", "#e575df",
    "#ef64a7"
)
names(mutation_pal) <- unique(unname(unlist(mapping_mutation_types_to_human_readable)))

ggsave_wrapper(
    show_palette(mutation_pal, "square", label_size = 4, font_family = "Arial"),
    plot_path("00_miscellaneous", "mutation_pal.svg"),
    width = 12, height = 4
)


# colors for KRAS alleles
allele_palette <- c(
    "KRAS G12A" = "#fb7810",
    "KRAS G12C" = "#0a4f4e",
    "KRAS G12D" = "#427ff5",
    "KRAS G12R" = "#af5fe4",
    "KRAS G12S" = "#f4cacb",
    "KRAS G12V" = "#f8d147",
    "KRAS G13C" = "#afe642",
    "KRAS G13D" = "#1c4585",
    "KRAS Q61H" = "#fd2c3b",
    "KRAS Q61K" = "#E76707",
    "KRAS Q61L" = "#56eead",
    "KRAS Q61R" = "#f33bea",
    "KRAS K117N" = "#E707C8",
    "KRAS A146T" = "#859947",
    "KRAS A146V" = "#9f75a7",
    "KRAS Other" = "grey75"
)

allele_palette <- c(allele_palette, "WT" = "grey50")
short_allele_pal <- allele_palette
names(short_allele_pal) <- str_remove_all(names(short_allele_pal), "KRAS ")

ggsave_wrapper(
    show_palette(short_allele_pal, "square", font_family = "Arial"),
    plot_path("00_miscellaneous", "short_alleles_pal.svg"),
    "small"
)


# colors for cancers
cancer_palette <- c(
    COAD = "darkolivegreen3",
    LUAD = "black",
    MM = "cornflowerblue",
    PAAD = "brown",
    SKCM = "darksalmon"
)

ggsave_wrapper(
    show_palette(cancer_palette, num_rows = 1, font_family = "Arial"),
    plot_path("00_miscellaneous", "cancer_palette.svg"),
    width = 6, height = 0.3
)



# colors for comutation vs mutual exclusivity
comut_mutex_pal <- c(
    comutation = "#1f71c2",  # blue
    exclusivity = "#3deb74"  # green
)
comut_updown_pal <- comut_mutex_pal
names(comut_updown_pal) <- c("increased", "reduced")


#' Convert from "comutation" and "exclusivity" to "increased" and "reduced".
switch_comut_terms <- function(old_terms) {
    conversion_vect <- names(comut_updown_pal)
    names(conversion_vect) <- names(comut_mutex_pal)
    return(unname(conversion_vect[old_terms]))
}


# Colors for extremes of synthetic lethal values.
synthetic_lethal_pal <- c(
    up = "#F47747",   # orange
    down = "#A33BD5"  # purple
)

ggsave_wrapper(
    show_palette(c(comut_mutex_pal, synthetic_lethal_pal),
                 font_family = "Arial"),
    plot_path("00_miscellaneous", "genetic_interaction_pal.svg"),
    width = 3.5, height = 2
)
