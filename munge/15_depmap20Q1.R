# Prepare the DepMap 2020Q1 data.

# Read in a DepMap data file.
# the `...` is passed to `read_csv()`.
read_depmap20Q1_csv <- function(x, ...) {
    file.path("data", "depmap20Q1", x) %>%
        read_csv(progress = FALSE, ...)
}



#### ---- Essential and nonessential genes ---- ####

cache("essentiality_tib",
{
    essentiality_tib <- bind_rows(
        read_depmap20Q1_csv("Achilles_common_essentials.csv") %>%
            add_column(label = "achilles_essential"),
        read_depmap20Q1_csv("common_essentials.csv") %>%
            add_column(label = "common_essential"),
        read_depmap20Q1_csv("nonessentials.csv") %>%
            add_column(label = "nonessential")
    ) %>%
        mutate(hugo_symbol = get_hugo_from_depmap_ids(gene),
               entrez_id = get_entrez_from_depmap_ids(gene)) %>%
        select(hugo_symbol, entrez_id, label)
    return(essentiality_tib)
})



#### ---- CNA ---- ####

cache("ccle_copy_number",
{
    ccle_copy_number <- read_depmap20Q1_csv("CCLE_gene_cn.csv") %>%
        rename(dep_map_id = X1) %>%
        pivot_longer(-dep_map_id,
                     names_to = "gene",
                     values_to = "log_copy_number") %>%
        mutate(hugo_symbol = get_hugo_from_depmap_ids(gene),
               entrez_id = get_entrez_from_depmap_ids(gene),
               copy_number = (2^log_copy_number) - 1,
               copy_number_label = case_when(copy_number < 0.5 ~ "homo_del",
                                             copy_number < 1.5 ~ "het_del",
                                             copy_number > 2.5 ~ "amp",
                                             TRUE ~ "norm")) %>%
        select(dep_map_id, hugo_symbol, entrez_id,
               log_copy_number, copy_number, copy_number_label)
    return(ccle_copy_number)
})



#### ---- Mutations ---- ####

ccle_mutations_cols <- cols(
    Hugo_Symbol = col_character(),
    Entrez_Gene_Id = col_character(),
    NCBI_Build = col_double(),
    Chromosome = col_character(),
    Start_position = col_double(),
    End_position = col_double(),
    Strand = col_character(),
    Variant_Classification = col_character(),
    Variant_Type = col_character(),
    Reference_Allele = col_character(),
    Tumor_Seq_Allele1 = col_character(),
    dbSNP_RS = col_character(),
    dbSNP_Val_Status = col_character(),
    Genome_Change = col_character(),
    Annotation_Transcript = col_character(),
    Tumor_Sample_Barcode = col_character(),
    cDNA_Change = col_character(),
    Codon_Change = col_character(),
    Protein_Change = col_character(),
    isDeleterious = col_logical(),
    isTCGAhotspot = col_logical(),
    TCGAhsCnt = col_double(),
    isCOSMIChotspot = col_logical(),
    COSMIChsCnt = col_double(),
    ExAC_AF = col_double(),
    CGA_WES_AC = col_character(),
    SangerWES_AC = col_character(),
    SangerRecalibWES_AC = col_character(),
    RNAseq_AC = col_character(),
    HC_AC = col_character(),
    RD_AC = col_character(),
    WGS_AC = col_character(),
    Variant_annotation = col_character(),
    DepMap_ID = col_character()
)

cache("ccle_mutations",
{
    ccle_mutations <- read_depmap20Q1_csv(
            "CCLE_mutations.csv", col_types = ccle_mutations_cols
        ) %>%
        janitor::clean_names() %>%
        select(
            dep_map_id,
            hugo_symbol, entrez_id = entrez_gene_id,
            ncbi_build:end_position,
            variant_annotation,
            variant_classification:tumor_seq_allele1,
            genome_change, codon_change:is_deleterious,
            is_tcga_hotspot = is_tcg_ahotspot,
            is_cosmic_hotspot = is_cosmi_chotspot
        )
    return(ccle_mutations)
})

cache("ccle_mutations_dmg", depends = "ccle_mutations",
{
    ccle_mutations_dmg <- ccle_mutations %>%
        filter(variant_annotation != "silent")
})



#### ---- KRAS mutations ---- ####

cache("ccle_kras_muts",
      depends = c("ccle_copy_number", "ccle_mutations_dmg"),
{
    ccle_kras_cna <- ccle_copy_number %>%
        filter(hugo_symbol == "KRAS" & copy_number_label == "amp") %>%
        select(dep_map_id, copy_number, copy_number_label)

    ccle_kras_muts <- ccle_mutations_dmg %>%
        filter(hugo_symbol == "KRAS") %>%
        select(dep_map_id, codon_change, protein_change) %>%
        mutate(
            aa_mod = str_remove(protein_change, "^p\\."),
            allele = ifelse(aa_mod %in% names(short_allele_pal),
                                 aa_mod, "other"),
            codon = ifelse(allele == "other", NA, str_extract(allele, "[:digit:]+"))
        ) %>%
        full_join(ccle_kras_cna, by = "dep_map_id") %>%
        mutate(
            codon_change = ifelse(is.na(codon_change), "WT", codon_change),
            protein_change = ifelse(is.na(protein_change), "WT", protein_change),
            aa_mod = ifelse(is.na(aa_mod), "WT", aa_mod),
            allele = ifelse(is.na(allele), "WT", allele),
            codon = ifelse(is.na(codon), "WT", codon),
            copy_number = ifelse(is.na(copy_number), 2, copy_number),
            copy_number_label = ifelse(is.na(copy_number_label),
                                       "norm", copy_number_label)
        )

    rm("ccle_kras_cna")
    return(ccle_kras_muts)
})



#### ---- RNA expression ---- ####

cache("ccle_expression",
{
    ccle_expression <- read_depmap20Q1_csv("CCLE_expression.csv") %>%
        rename(dep_map_id = X1) %>%
        pivot_longer(-dep_map_id,
                     names_to = "gene",
                     values_to = "rna_expression") %>%
        mutate(hugo_symbol = get_hugo_from_depmap_ids(gene),
               entrez_id = get_entrez_from_depmap_ids(gene)) %>%
        select(dep_map_id, hugo_symbol, entrez_id, rna_expression)
    return(ccle_expression)
})



#### ---- Sample (cell line) Information ---- ####

ccle_cell_lines_cols <- cols(
    DepMap_ID = col_character(),
    stripped_cell_line_name = col_character(),
    CCLE_Name = col_character(),
    alias = col_character(),
    COSMIC_ID = col_double(),
    lineage = col_character(),
    lineage_subtype = col_character(),
    lineage_sub_subtype = col_character(),
    lineage_molecular_subtype = col_character(),
    sex = col_character(),
    source = col_character(),
    Achilles_n_replicates = col_double(),
    cell_line_NNMD = col_double(),
    culture_type = col_character(),
    culture_medium = col_character(),
    cas9_activity = col_character(),
    RRID = col_character(),
    sample_collection_site = col_character(),
    primary_or_metastasis = col_character(),
    disease = col_character(),
    disease_subtype = col_character(),
    age = col_double(),
    Sanger_model_ID = col_character(),
    additional_info = col_character()
)

cache("ccle_cell_lines",
{
    ccle_cell_lines <- read_depmap20Q1_csv(
            "sample_info_v2.csv", col_types = ccle_cell_lines_cols
        ) %>%
        janitor::clean_names() %>%
        mutate_at(c("lineage", "lineage_subtype",
                    "lineage_sub_subtype", "lineage_molecular_subtype",
                    "disease", "disease_subtype"),
                  str_to_lower) %>%
        mutate(
            cancer = case_when(
                lineage_subtype == "colorectal_adenocarcinoma" ~ "COAD",
                str_detect(disease_subtype, "carcinoma") &
                    disease == "colon/colorectal cancer" ~ "COAD",
                disease == "pancreatic cancer" &
                    str_detect(disease_subtype, "adenocarcinoma") ~ "PAAD",
                disease == "pancreatic cancer" &
                    disease_subtype == "carcinoma" ~ "PAAD",
                disease == "myeloma" ~ "MM",
                lineage_sub_subtype == "nsclc_adenocarcinoma" ~ "LUAD",
                disease == "lung cancer" &
                    str_detect(disease_subtype, "nsclc") &
                    str_detect(disease_subtype,
                               "adenocarcinoma|unspecified") ~ "LUAD",
                TRUE ~ NA_character_
        ))
    return(ccle_cell_lines)
})



#### ---- Gene effect ---- ####

cache("gene_effect",
{
    gene_effect <- read_depmap20Q1_csv("Achilles_gene_effect.csv") %>%
        rename(dep_map_id = X1) %>%
        pivot_longer(-dep_map_id,
                     names_to = "gene",
                     values_to = "gene_effect") %>%
        mutate(hugo_symbol = get_hugo_from_depmap_ids(gene),
               entrez_id = get_entrez_from_depmap_ids(gene)) %>%
        select(dep_map_id, hugo_symbol, entrez_id, gene_effect)
    return(gene_effect)
})
