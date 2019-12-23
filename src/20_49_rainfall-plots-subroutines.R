# Data frames and subroutines for preparing MAFs from cancer data.
# Accompanies "src/20_50_rainfall-plots.R"


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


# Put Variant_Calssification values to match those expected for MAF.
fix_variant_classification <- function(vc) {
    new_vc <- str_replace_us(vc) %>%
        str_to_title() %>%
        str_replace_sp()

    new_vc[new_vc == "Inframe_Deletion"] <- "In_Frame_Del"
    new_vc[new_vc == "Inframe_Insertion"] <- "In_Frame_Ins"
    new_vc[new_vc == "Frameshift_Deletion"] <- "Frame_Shift_Del"
    new_vc[new_vc == "Frameshift_Insertion"] <- "Frame_Shift_Ins"

    return(new_vc)
}

# Add the MAF-required columns.
cancer_mut_df_to_maf <- function(df) {
    df %>%
        group_by(tumor_sample_barcode) %>%
        dplyr::mutate(
            Hugo_Symbol = hugo_symbol,
            Chromosome = extract_chromosome(genomic_position),
            Start_Position = extract_genomic_position(genomic_position),
            End_Position = Start_Position,
            Reference_Allele = extract_reference_allele(genomic_position),
            Tumor_Seq_Allele2 = extract_tumor_allele(genomic_position),
            Variant_Classification = fix_variant_classification(mutation_type),
            Variant_Type = "SNP",
            Tumor_Sample_Barcode = tumor_sample_barcode,
            Protein_Change = paste0("p.", amino_acid_change)
        ) %>%
        ungroup()
}




# A tibble of cancer mutation data that is MAF-compliant.
ProjectTemplate::cache("cancer_coding_muts_maf", depends = "cancer_muts_df",
{
    cancer_coding_muts_maf <- cancer_coding_muts_df %>%
        filter(!is_hypermutant & cancer != "SKCM") %>%
        cancer_mut_df_to_maf() %>%
        filter(!is.na(Hugo_Symbol))
    return(cancer_coding_muts_maf)
})


# A tibble of cancer mutation data that is MAF-compliant.
# DO NOT USE THIS FOR THE RAINFALL PLOTS! -- only lollipop plots
ProjectTemplate::cache("cancer_full_coding_muts_maf",
                       depends = "cancer_full_coding_muts_df",
{
    cancer_full_coding_muts_maf <- cancer_full_coding_muts_df %>%
        filter(!is_hypermutant & cancer != "SKCM") %>%
        cancer_mut_df_to_maf() %>%
        filter(!is.na(Hugo_Symbol))
    return(cancer_full_coding_muts_maf)
})



# Filter the complete MAF data frame and return a MAF object for 'maftools'.
# If `replace_kras_with_allele`: set Hugo_Symbol as the allele for KRAS.
# If `group_other_alleles`: then all other KRAS alleles become "KRAS_other".
get_maf <- function(for_genes, in_cancer, kras_allele,
                    replace_kras_with_allele = FALSE,
                    group_other_alleles = TRUE) {
    df <- cancer_coding_muts_maf %>%
        dplyr::filter(cancer %in% !!in_cancer & hugo_symbol %in% !!for_genes)

    if (replace_kras_with_allele) {

        if (group_other_alleles) {
            df <- df %>%
                mutate(ras_allele = ifelse(
                    ras_allele == !!kras_allele, ras_allele, "KRAS_other"
                ))
        }

        df <- df %>%
            mutate(Hugo_Symbol = ifelse(
                Hugo_Symbol == "KRAS", str_replace_us(ras_allele), Hugo_Symbol
            ))
    }

    maf <- df %>%
        dplyr::select(
            Hugo_Symbol, Chromosome, Start_Position, End_Position,
            Reference_Allele, Tumor_Seq_Allele2, Variant_Classification,
            Variant_Type, Tumor_Sample_Barcode, VAF, Protein_Change
        ) %>%
        as.data.frame()

    clinical_data <- cancer_coding_muts_maf %>%
        dplyr::filter(cancer %in% !!in_cancer & hugo_symbol %in% !!for_genes) %>%
        select(Tumor_Sample_Barcode, ras_allele) %>%
        dplyr::rename(kras_allele = "ras_allele") %>%
        mutate(
            kras_allele = str_remove(kras_allele, "KRAS_"),
            kras_allele = fct_lump(kras_allele, prop = 0.01)
        ) %>%
        unique() %>%
        as.data.frame()
    read.maf(maf, clinicalData = clinical_data, verbose = FALSE)
}

