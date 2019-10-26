# Rainfall plots for the more interesting and strongest genetic interactions.

# CAREFUL: THIS BRINGS IN A LOT OF BIOCONDUCTOR BAGGAGE!
library(maftools)


#### ---- Prepare data for 'maftools' ---- ####
# It seems worth-while to prepare data for use with the 'maftools' pacakge.
# As the name implies, the input is usually data in the Mutation Annotation
# Format (MAF).
#
# --- MAF ---
# Required fields:
#   Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele,
#   Tumor_Seq_Allele2, Variant_Classification, Variant_Type,
#   Tumor_Sample_Barcode.
#
# Recommended fields:
#   VAF, Protein_Change

# Get chromosome.
extract_chromosome <- function(genomic_position) {
    str_split_fixed(genomic_position, ":", 2)[, 1] %>%
        unlist()
}

# Get position on the chromosome.
extract_genomic_position <- function(genomic_position) {
    str_split_fixed(genomic_position, ":", 3)[, 2] %>%
        unlist() %>%
        as.numeric()
}

# Get reference allele.
extract_reference_allele <- function(genomic_position) {
    acgt <- str_split_fixed(genomic_position, ":", 3)[, 3] %>% unlist()
    str_split_fixed(acgt, ",", 2)[, 1]
}

# Get mutation allele.
extract_tumor_allele <- function(genomic_position) {
    acgt <- str_split_fixed(genomic_position, ":", 3)[, 3] %>% unlist()
    str_split_fixed(acgt, ",", 2)[, 2]
}


ProjectTemplate::cache("cancer_muts_maf", depends = "cancer_muts_df",
{
    cancer_muts_maf <- cancer_muts_df %>%
        dplyr::mutate(
            Hugo_Symbol = hugo_symbol,
            Chromosome = extract_chromosome(genomic_position),
            Start_Position = extract_genomic_position(genomic_position),
            End_Position = Start_Position,
            Reference_Allele = extract_reference_allele(genomic_position),
            Tumor_Seq_Allele2 = extract_tumor_allele(genomic_position),
            Variant_Classification = mutation_type,
            Variant_Type = mutation_type,
            Tumor_Sample_Barcode = tumor_sample_barcode,
            Protein_Change = paste0("p.", amino_acid_change)
        )
    return(cancer_muts_maf)
})


# Filter the complete MAF data frame and return a MAF object for 'maftools'
get_maf <- function(for_genes, in_cancer) {
    browser()
    cancer_muts_maf %>%
        dplyr::filter(cancer %in% !!in_cancer & hugo_symbol %in% !!for_genes) %>%
        dplyr::select(
            Hugo_Symbol, Chromosome, Start_Position, End_Position,
            Reference_Allele, Tumor_Seq_Allele2, Variant_Classification,
            Variant_Type, Tumor_Sample_Barcode, VAF, Protein_Change
        ) %>%
        dplyr::as._data_frame() %>%
        dplyr::read.maf()
}



