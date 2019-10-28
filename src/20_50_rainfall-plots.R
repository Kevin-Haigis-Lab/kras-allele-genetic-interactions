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


ProjectTemplate::cache("cancer_coding_muts_maf", depends = "cancer_muts_df",
{
    cancer_coding_muts_maf <- cancer_coding_muts_df %>%
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
        ungroup() %>%
        filter(!is.na(Hugo_Symbol))
    return(cancer_coding_muts_maf)
})


# Filter the complete MAF data frame and return a MAF object for 'maftools'
# If `replace_kras_with_allele`: set Hugo_Symbol as the allele for KRAS
# If `group_other_alleles`: then all other KRAS alleles become "KRAS_other"
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




#### ---- Oncoplots for interactors of KRAS by allele ---- ####
# Plot the `maftools::oncostrip()` for the genes mutually exclusive with one
# allele for all of the alleles.


# Get all KRAS alleles for a cancer
get_alleles_in_cancer <- function(cancer) {
    cancer_coding_muts_df %>%
        filter(cancer == !!cancer) %>%
        u_pull(ras_allele)
}


# Sort genes into: KRAS, "KRAS other", by mut. freq.
sort_genes <- function(maf) {
    gene_order <- maf@data %>%
        as_tibble() %>%
        count(Hugo_Symbol) %>%
        arrange(-n) %>%
        pull(Hugo_Symbol)

    kras_allele <- gene_order[str_detect(gene_order, "KRAS ")]
    if (any(kras_allele == "KRAS other")) {
        kras_allele <- c(
            kras_allele[kras_allele != "KRAS other"],
            "KRAS other"
        )
    }
    gene_order <- c(
        kras_allele,
        gene_order[!(gene_order %in% kras_allele)]
    )
    return(gene_order)
}


# If there is a KRAS allele in `genes` return it, otherwise return "WT".
any_kras_allele <- function(genes) {
    kras_gene <- genes[str_detect(genes, "KRAS ")]
    if (length(kras_gene) > 0) { return(kras_gene[[1]]) }
    else { return("WT") }
}


# Sort samples using a weighted metric.
#  The first weight is by KRAS allele, with "KRAS_other" at the end.
#  The second is by the order of the genes in `ordered_genes`.
sort_samples <- function(maf, ordered_genes) {
    gene_score_tib <- tibble(Hugo_Symbol = rev(ordered_genes)) %>%
        mutate(value = 2 ^ ((1:n()) - 1))

    kras_tib <- tibble(
        kras_allele = rev(ordered_genes[str_detect(ordered_genes, "KRAS ")])
    ) %>%
        mutate(kras_value = 1:n())

    maf@data %>%
        as_tibble() %>%
        select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
        left_join(gene_score_tib, by = "Hugo_Symbol") %>%
        group_by(Tumor_Sample_Barcode) %>%
        summarise(score = sum(value),
                  kras_allele = any_kras_allele(Hugo_Symbol)) %>%
        ungroup() %>%
        left_join(kras_tib, by = "kras_allele") %>%
        mutate(kras_value = ifelse(is.na(kras_value), 0, kras_value)) %>%
        arrange(-kras_value, -score) %>%
        pull(Tumor_Sample_Barcode)
}


# Sort the genes for oncostrip and keep KRAS at the top.
sort_maf <- function(maf) {
    gene_order <- sort_genes(maf)
    sample_order <- sort_samples(maf, ordered_genes = gene_order)
    return(list(
        gene_order = gene_order,
        sample_order = sample_order
    ))
}


# MAIN PLOTTING FUNCTION
# Plot the oncoplot for a `cancer` and KRAS `allele` for the
#   top `top_n_genes` genes, ignoring any in `ignore_genes`.
# (The `...` doesn't go anywhere.)
allele_exclusivity_oncostrip <- function(cancer, allele,
                                         save_name_template,
                                         interaction_type = "exclusivity",
                                         top_n_genes = 20,
                                         ignore_genes = c(),
                                         annotate_kras_allele = FALSE,
                                         ...) {
    mutex_genes <- genetic_interaction_gr %N>%
        filter(!name %in% !!ignore_genes) %E>%
        filter(cancer == !!cancer &
               kras_allele == !!allele &
               genetic_interaction %in% !!interaction_type)  %>%
        top_n(top_n_genes, -p_val) %N>%
        filter(centrality_degree(mode = "all") > 0) %N>%
        as_tibble() %>%
        filter(!is_kras) %>%
        pull(name)

    if (length(mutex_genes) == 0) return()

    mutex_genes <- c(
        mutex_genes,
        "KRAS"
    )

    maf <- get_maf(mutex_genes, cancer, allele,
                   replace_kras_with_allele = TRUE,
                   group_other_alleles = TRUE)

    gene_and_sample_orders <- sort_maf(maf)

    allele_short <- str_remove(allele, "KRAS_")
    save_path <- plot_path(
        "20_50_rainfall-plots",
        glue(save_name_template)
    )

    if (annotate_kras_allele) {
        pal <- short_allele_pal[names(short_allele_pal) %in% getClinicalData(maf)$kras_allele]

        svg(save_path, width = 9, height = 5)
        oncoplot(
            maf,
            # genes = gene_and_sample_orders$gene_order[!str_detect(gene_and_sample_orders$gene_order, "KRAS ")],
            genes = gene_and_sample_orders$gene_order,
            sampleOrder = gene_and_sample_orders$sample_order,
            top = top_n_genes,
            keepGeneOrder = TRUE,
            GeneOrderSort = FALSE,
            removeNonMutated = FALSE,
            clinicalFeatures = 'kras_allele',
            annotationColor = list(kras_allele = pal),
            sepwd_samples = 0,
            fontSize = 0.6
        )
        dev.off()
    } else {
        svg(save_path, width = 9, height = 5)
        oncoplot(
            maf,
            genes = gene_and_sample_orders$gene_order,
            sampleOrder = gene_and_sample_orders$sample_order,
            top = top_n_genes,
            keepGeneOrder = TRUE,
            GeneOrderSort = FALSE,
            sepwd_samples = 0,
            fontSize = 0.6
        )
        dev.off()
    }
}


# Uninteresting genes
uninteresting_genes <- c("TTN", "FN1")

# Genes to ignore in most situations.
commonly_interacting_genes <- c(
    uninteresting_genes,
    "BRAF", "APC", "TP53", "NRAS", "PIK3CA"
)


cancer_counts_df <- cancer_muts_df %>%
    group_by(cancer, ras_allele, tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    count(cancer, ras_allele)

cancer_counts_df %>%
    filter(n > 10) %>%
    dplyr::rename(allele = "ras_allele") %>%
    pwalk(
        allele_exclusivity_oncostrip,
        save_name_template = "{cancer}_{allele_short}_exclusivity_oncostrip_notCommonInteractors.svg",
        top_n_genes = 15,
        ignore_genes = commonly_interacting_genes,
        annotate_kras_allele = TRUE
    )


cancer_counts_df %>%
    filter(n > 10) %>%
    dplyr::rename(allele = "ras_allele") %>%
    pwalk(
        allele_exclusivity_oncostrip,
        save_name_template = "{cancer}_{allele_short}_exclusivity_oncostrip_allInteractors.svg",
        top_n_genes = 15,
        ignore_genes = uninteresting_genes,
        annotate_kras_allele = TRUE
    )

cancer_counts_df %>%
    filter(n > 10) %>%
    dplyr::rename(allele = "ras_allele") %>%
    pwalk(
        allele_exclusivity_oncostrip,
        save_name_template = "{cancer}_{allele_short}_bothInteractions_oncostrip_allInteractors.svg",
        interaction_type = c("comutation", "exclusivity"),
        top_n_genes = 15,
        ignore_genes = uninteresting_genes,
        annotate_kras_allele = TRUE
    )

cancer_counts_df %>%
    filter(n > 10) %>%
    dplyr::rename(allele = "ras_allele") %>%
    pwalk(
        allele_exclusivity_oncostrip,
        save_name_template = "{cancer}_{allele_short}_comutation_oncostrip_allInteractors.svg",
        interaction_type = "comutation",
        top_n_genes = 15,
        ignore_genes = uninteresting_genes,
        annotate_kras_allele = TRUE
    )
