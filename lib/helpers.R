
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

# helpful data for identifying KRAS mutations
kras_hotspot_codons <- list(num = c(12, 13, 61, 117, 146))
kras_hotspot_codons$char <- as.character(kras_hotspot_codons$num)
kras_hotspot_codons$regex <- paste0(kras_hotspot_codons$char, collapse = "|")

# coding mutation variant classification terms
coding_mut_var_classes <- c(
    "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
    "Nonstop_Mutation", "Splice_Site", "Start_Codon_Del", "Start_Codon_Ins",
    "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins"
)

