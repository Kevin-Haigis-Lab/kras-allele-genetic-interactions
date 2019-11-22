# Rainfall plots for the more interesting and strongest genetic interactions.


#### ---- Rate of mutation for each gene per allele, WT, and overall ---- ####

ProjectTemplate::cache("comutation_rates_df",
                       depends = "cancer_coding_muts_df",
{
    cancer_mut_comutation_df <- cancer_coding_muts_df %>%
        filter(!is_hypermutant) %>%
        group_by(cancer) %>%
        mutate(number_samples_cancer = n_distinct(tumor_sample_barcode)) %>%
        group_by(cancer, hugo_symbol, number_samples_cancer) %>%
        summarise(number_mutations_cancer = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        mutate(freq_mutations_cancer = number_mutations_cancer / number_samples_cancer)

    kras_mut_comutation_df <- cancer_coding_muts_df %>%
        filter(!is_hypermutant) %>%
        group_by(cancer, ras) %>%
        mutate(number_samples_ras = n_distinct(tumor_sample_barcode)) %>%
        group_by(cancer, ras, hugo_symbol, number_samples_ras) %>%
        summarise(number_mutations_ras = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        mutate(freq_mutations_ras = number_mutations_ras / number_samples_ras) %>%
        dplyr::rename(ras_allele = "ras")

    allele_mut_comutation <- cancer_coding_muts_df %>%
        filter(!is_hypermutant) %>%
        group_by(cancer, ras_allele) %>%
        mutate(number_samples_ras = n_distinct(tumor_sample_barcode)) %>%
        group_by(cancer, ras_allele, hugo_symbol, number_samples_ras) %>%
        summarise(number_mutations_ras = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        mutate(freq_mutations_ras = number_mutations_ras / number_samples_ras)


    comutation_rates_df <- full_join(
            cancer_mut_comutation_df,
            bind_rows(kras_mut_comutation_df, allele_mut_comutation),
            by = c("cancer", "hugo_symbol")
        ) %>%
        group_by(hugo_symbol) %>%
        unique() %>%
        ungroup()

    return(comutation_rates_df)
})


# comutation_rates_df %>%
#     filter(cancer == "PAAD" & hugo_symbol == "TP53") %>%
#     filter(ras_allele %in% c("KRAS_G12D", "KRAS_G12V", "KRAS_Q61H", "WT", "KRAS")) %>%
#     select(ras_allele, freq_mutations_cancer, freq_mutations_ras)


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
        kras_allele = rev(ordered_genes[str_detect(ordered_genes, "KRAS")])
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


# Tibble of number of samples per cancer to be used for the oncoplot.
cancer_count_tib <- cancer_muts_df %>%
    group_by(cancer, tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    unique() %>%
    count(cancer)


oncoplot_wrapper <- function(maf,
                             gene_and_sample_orders,
                             save_path,
                             top = 10,
                             cancer = NULL,
                             annotate_kras_allele = FALSE,
                             save_size = list(width = 9, height = 5)) {

    if (!is.null(cancer)) {
        cohortSize <- unlist(filter(cancer_count_tib, cancer == !!cancer)$n)
    } else {
        cohortSize <- NULL
    }

    if (annotate_kras_allele) {
        idx <- names(short_allele_pal) %in% getClinicalData(maf)$kras_allele
        pal <- short_allele_pal[idx]

        svg(save_path, width = save_size$width, height = save_size$height)
        oncoplot(
            maf,
            genes = gene_and_sample_orders$gene_order,
            sampleOrder = gene_and_sample_orders$sample_order,
            top = top,
            cohortSize = cohortSize,
            keepGeneOrder = TRUE,
            GeneOrderSort = FALSE,
            removeNonMutated = TRUE,
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
            top = top,
            cohortSize = cohortSize,
            keepGeneOrder = TRUE,
            GeneOrderSort = FALSE,
            removeNonMutated = TRUE,
            sepwd_samples = 0,
            fontSize = 0.6
        )
        dev.off()
    }
}


# MAIN PLOTTING FUNCTION
# Plot the oncoplot for a `cancer` and KRAS `allele` for the
#   top `top_n_genes` genes, ignoring any in `ignore_genes`.
# (The `...` doesn't go anywhere.)
allele_exclusivity_oncoplot <- function(cancer, allele,
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

    oncoplot_wrapper(
        maf = maf,
        gene_and_sample_orders = gene_and_sample_orders,
        save_path = save_path,
        top = top_n_genes,
        cancer = cancer,
        annotate_kras_allele = annotate_kras_allele
    )
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
    count(cancer, ras_allele) %>%
    dplyr::rename(allele = "ras_allele")

cancer_counts_df %>%
    filter(n > 10) %>%
    pwalk(
        allele_exclusivity_oncoplot,
        save_name_template = "{cancer}_{allele_short}_exclusivity_oncostrip_notCommonInteractors.svg",
        top_n_genes = 15,
        ignore_genes = commonly_interacting_genes,
        annotate_kras_allele = TRUE
    )

cancer_counts_df %>%
    filter(n > 10) %>%
    pwalk(
        allele_exclusivity_oncoplot,
        save_name_template = "{cancer}_{allele_short}_exclusivity_oncostrip_allInteractors.svg",
        top_n_genes = 15,
        ignore_genes = uninteresting_genes,
        annotate_kras_allele = TRUE
    )

cancer_counts_df %>%
    filter(n > 10) %>%
    pwalk(
        allele_exclusivity_oncoplot,
        save_name_template = "{cancer}_{allele_short}_bothInteractions_oncostrip_allInteractors.svg",
        interaction_type = c("comutation", "exclusivity"),
        top_n_genes = 15,
        ignore_genes = uninteresting_genes,
        annotate_kras_allele = TRUE
    )

cancer_counts_df %>%
    filter(n > 10) %>%
    pwalk(
        allele_exclusivity_oncoplot,
        save_name_template = "{cancer}_{allele_short}_comutation_oncostrip_allInteractors.svg",
        interaction_type = "comutation",
        top_n_genes = 15,
        ignore_genes = uninteresting_genes,
        annotate_kras_allele = TRUE
    )



#### ---- Specific genes to  make oncoplots for ---- ####



# MAIN PLOTTING FUNCTION
# Plot the oncoplot for a `cancer` and KRAS `allele` for a specific set of
#   genes. The genes are checked against the interaction graph.
# (The `...` doesn't go anywhere.)
allele_specificgenes_oncoplot <- function(cancer,
                                          kras_allele,
                                          genes,
                                          interaction_type,
                                          save_name,
                                          save_dir,
                                          ignore_genes = c(),
                                          top_n_genes = 10,
                                          replace_kras_with_allele = TRUE,
                                          annotate_kras_allele = FALSE,
                                          ...) {
    genes <- unlist(genes)

    checked_genes <- genetic_interaction_gr %N>%
        filter(!(name %in% !!ignore_genes)) %E>%
        filter(cancer == !!cancer &
               (kras_allele == !!kras_allele | !!kras_allele == "KRAS_ALL")) %N>%
        filter(centrality_degree(mode = "all") > 0) %>%
        as_tibble() %>%
        filter(!is_kras & name %in% !!genes) %>%
        u_pull(name)

    cat(glue("{cancer}, {kras_allele}: {n_distinct(genes) - n_distinct(checked_genes)} gene(s) were removed."), "\n")

    if (length(checked_genes) == 0) {
        cat("**There were zero genes to plot; returning early.\n")
        return()
    }

    checked_genes <- c(
        checked_genes,
        "KRAS"
    )

    maf <- get_maf(checked_genes, cancer, kras_allele,
                   replace_kras_with_allele = replace_kras_with_allele,
                   group_other_alleles = replace_kras_with_allele)

    gene_and_sample_orders <- sort_maf(maf)

    save_path <- plot_path(
        save_dir,
        save_name
    )

    img_height <- ((length(checked_genes) + 1) / 5) + 2

    oncoplot_wrapper(
        maf = maf,
        gene_and_sample_orders = gene_and_sample_orders,
        save_path = save_path,
        top = top_n_genes + 1 + as.numeric(replace_kras_with_allele),
        cancer = cancer,
        annotate_kras_allele = annotate_kras_allele,
        save_size = list(width = 9, height = img_height)
    )
}


specific_oncoplot_info_tib <- bind_rows(
    tibble(
        cancer = "COAD",
        allele = "A146T",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("PORCN", "RGS20", "DAPK1"),
    ),
    tibble(
        cancer = "COAD",
        allele = "G12D",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("SCN10A", "DICER1", "SDK1", "TRPM2")
    ),
    tibble(
        cancer = "COAD",
        allele = "G12D",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("MAGEC1", "AMER1", "TGIF1")
    ),
    tibble(
        cancer = "COAD",
        allele = "G12V",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("DNAH2", "RYR1", "PKHD1", "DIDO1", "PKD1", "KALRN")
    ),
    tibble(
        cancer = "COAD",
        allele = "G12V",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("APC", "PIK3CA", "SMAD4", "AMER1", "MCC")
    ),
    tibble(
        cancer = "COAD",
        allele = "G13D",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("ANK3", "KIF4B", "SPHKAP")
    ),
    tibble(
        cancer = "COAD",
        allele = "G13D",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("TP53", "CAMTA1", "ZFHX4")
    ),
    tibble(
        cancer = "LUAD",
        allele = "G12C",
        interaction_type = "exclusivity",
        name_suffix = "(MUC)",
        genes = c("MUC4", "MUC17")
    ),
    tibble(
        cancer = "LUAD",
        allele = "G12C",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("CSMD2", "DNAH5", "SMARCA4", "CHD5"),
    ),
    tibble(
        cancer = "LUAD",
        allele = "G12C",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("PRDM9", "RIMS2", "STK11", "SLITRK2")
    ),
    tibble(
        cancer = "LUAD",
        allele = "G12D",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("TP53")
    ),
    tibble(
        cancer = "LUAD",
        allele = "G12V",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("FAT4", "LAMA2", "SPHKAP", "ZNF804A")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12D",
        interaction_type = "exclusivity",
        name_suffix = "(RYR)",
        genes = c("RYR3", "RYR2")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12D",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("RYR3", "RYR2", "TGFBR2", "GNAS")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12D",
        interaction_type = "comutation",
        name_suffix = "_TP53",
        genes = c("TP53")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12D",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("ABCC9", "ZNF831")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12R",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("RNF43", "DNAH11", "GNAS", "DNAH5", "ARID2", "SMARCA4", "TAF1")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12R",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("RAF1", "IRAK1")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12V",
        interaction_type = "comutation",
        name_suffix = "_TP53",
        genes = c("TP53")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12V",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("RYR1", "FAT4", "DNAH17", "MYH1", "LRRK2")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12V",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("RYR3", "PCDHA4", "MAGI1")
    ),
    tibble(
        cancer = "PAAD",
        allele = "G12V",
        interaction_type = "exclusivitycomutation",
        name_suffix = "(RYR)",
        genes = c("RYR1", "RYR3")
    ),
    tibble(
        cancer = "PAAD",
        allele = "Q61H",
        interaction_type = "exclusivity",
        name_suffix = "",
        genes = c("SMAD4", "ARID1A", "DNAH5", "ATM", "KDM6A")
    ),
    tibble(
        cancer = "PAAD",
        allele = "Q61H",
        interaction_type = "comutation",
        name_suffix = "",
        genes = c("TP53", "PTPRB", "CREB3L3")
    ),
    tibble(
        cancer = "PAAD",
        allele = "Q61H",
        interaction_type = "comutation",
        name_suffix = "_TP53",
        genes = c("TP53")
    ),
    tibble(
        cancer = "PAAD",
        allele = "Q61H",
        interaction_type = "exclusivitycomutation",
        name_suffix = "_oncogenes",
        genes = c("TP53", "SMAD4", "ARID1A")
    ),
    # tibble(
    #     cancer = "",
    #     allele = "",
    #     interaction_type = "exclusivity comutation",
    #     name_suffix = "",
    #     genes = c()
    # )
) %>%
    group_by(cancer, allele, interaction_type, name_suffix) %>%
    summarise(genes = list(genes)) %>%
    ungroup() %>%
    mutate(kras_allele = paste0("KRAS_", allele))

specific_oncoplot_info_tib %>%
    mutate(
        save_name = paste0(cancer, "_",
                           allele, "_",
                           interaction_type,
                           "_oncostrip_select",
                           name_suffix,
                           ".svg")
    ) %>%
    pwalk(
        allele_specificgenes_oncoplot,
        top_n_genes = 15,
        save_dir = "20_50_rainfall-plots-select",
        annotate_kras_allele = TRUE
    )


#### ---- Oncoplots for known interactors ---- ####

known_oncoplot_info_tib <- bind_rows(
    tibble(
        cancer = "COAD",
        interaction_type = "comutation",
        name_suffix = "_KNOWN",
        allele = "ALL",
        genes = c("PIK3CA")
    ),
    tibble(
        cancer = "COAD",
        interaction_type = "exclusivity",
        name_suffix = "_KNOWN",
        allele = "ALL",
        genes = c("BRAF", "NRAS")
    ),
    tibble(
        cancer = "LUAD",
        interaction_type = "exclusivity",
        name_suffix = "_KNOWN",
        allele = "ALL",
        genes = c("BRAF", "EGFR")
    ),
    tibble(
        cancer = "LUAD",
        interaction_type = "comutation",
        name_suffix = "_KNOWN",
        allele = "ALL",
        genes = c("STK11")
    ),
    tibble(
        cancer = "MM",
        interaction_type = "exclusivity",
        name_suffix = "_KNOWN",
        allele = "ALL",
        genes = c("NRAS", "BRAF")
    ),
    tibble(
        cancer = "PAAD",
        interaction_type = "exclusivitycomutation",
        name_suffix = "_KNOWN",
        allele = "ALL",
        genes = c("BRAF", "TP53")
    )
) %>%
    group_by(cancer, allele, interaction_type, name_suffix) %>%
    summarise(genes = list(genes)) %>%
    ungroup() %>%
    mutate(kras_allele = paste0("KRAS_", allele))

known_oncoplot_info_tib %>%
    mutate(
        save_name = paste0(cancer, "_",
                           allele, "_",
                           interaction_type,
                           "_oncostrip_select",
                           name_suffix,
                           ".svg")
    ) %>%
    pwalk(
        allele_specificgenes_oncoplot,
        top_n_genes = 15,
        save_dir = "20_50_rainfall-plots-select",
        replace_kras_with_allele = FALSE,
        annotate_kras_allele = TRUE
    )


#### ---- Oncoplots for a priori lists ---- ####
# Oncoplots to accompany the subsetted genetic networks.

gene_sets <- list(
    kegg = unique(kegg_geneset_df$hugo_symbol),
    cgc = unique(filter(cosmic_cgc_df, tier == 1)$hugo_symbol),
    bioid = unique(kras_interactors_bioid_df$hugo_symbol)
)

for (i in 1:length(gene_sets)) {
    gene_set_name <- names(gene_sets)[[i]]
    gene_set <- unname(unlist(gene_sets[i]))


    cancer_counts_df %>%
        mutate(kras_allele = allele,
               allele = str_remove_all(allele, "KRAS_"),
               save_name = paste0(cancer, "_",
                                  allele, "_",
                                  !!gene_set_name, "_",
                                  "oncoplot.svg")) %>%
        pwalk(allele_specificgenes_oncoplot,
              genes = gene_set,
              top_n_genes = 10,
              save_dir = "20_50_rainfall-plots-apriori",
              annotate_kras_allele = TRUE)
}
