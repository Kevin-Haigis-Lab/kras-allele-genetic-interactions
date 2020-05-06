
# Prepare the data to be used to model depletion effects from the KRAS allele


################################################################################
################################################################################

# NOTES: this data is old and this script is for legacy purposes, only.
# The munge script "16_depmap20Q1_model-data-prep.R" uses the newest available
# data.

# To run this script, set `RUN_LEGACY_DEPMAP_MUNGE = TRUE`
# Everything is wrapped in a massive if-statement.

################################################################################
################################################################################

RUN_LEGACY_DEPMAP_MUNGE <- FALSE
if (RUN_LEGACY_DEPMAP_MUNGE) {


#### ---- General Data ---- ####

# Cell lines with two oncogenic KRAS mutations
info(logger, "Preparing a vector of KRAS double mutants.")
kras_double_mutants <- kras_mutation_tib %>%
    filter(is_hotspot) %>%
    group_by(dep_map_id) %>%
    filter(n_distinct(allele) > 1) %>%
    ungroup() %>%
    pull(dep_map_id) %>%
    unique()
info(logger, glue("There were {n_distinct(kras_double_mutants)} KRAS double mutants."))
cache("kras_double_mutants")


info(logger, "Preparing a vector of BRAF, NRAS, and EGFR mutants.")
mapk_muts <- ccle_mutations_coding %>%
    filter(
        (hugo_symbol == "NRAS" & str_detect(protein_change, "G12|G13|Q61")) |
        (hugo_symbol == "BRAF" & str_detect(protein_change, "V600")) |
        (hugo_symbol == "EGFR" & is_cosmic_hotspot)
    ) %>%
    pull(dep_map_id) %>%
    unique()
info(logger, glue("There were {n_distinct(mapk_muts)} MAPK mutants."))
cache("mapk_muts")



# a tibble of oncogenic KRAS mutations to merge with modeling data
kras_mut_tib_tomerge <- kras_mutation_tib %>%
    filter(
        is_hotspot & !(dep_map_id %in% c(kras_double_mutants, mapk_muts))
    ) %>%
    select(dep_map_id, codon, allele, copy_number_label) %>%
    unique()

# check that each sample only has one KRAS mutation to merge with gene effect
n_kras_muts <- kras_mut_tib_tomerge %>%
    group_by(dep_map_id) %>%
    filter(n() > 1) %>%
    pull(dep_map_id)

if (length(n_kras_muts) == 0) {
    info(logger, "All samples have 0 or 1 KRAS mutations being added to the data.")
} else {
    error(logger, "There are cell lines with multiple KRAS mutations being merged.")
}


# Get essential and nonessential genes.
essential_genes <- essentiality_tib %>%
    filter(label == "achilles_essential") %>%
    jhcutils::u_pull(hugo_symbol)

nonessential_genes <- essentiality_tib %>%
    filter(label == "nonessential") %>%
    jhcutils::u_pull(hugo_symbol)

sample_info_tomerge <- cell_lines %>%
    select(dep_map_id, cancer) %>%
    filter(!is.na(cancer)) %>%
    unique()

stopifnot(!any(table(sample_info_tomerge$dep_map_id) > 1))



# Combination of CNA and coding mutation alterations.
info(logger, "Collection CNA and coding mutations for each gene and cell line.")
gene_is_aletered_tib <- bind_rows(
    {
        ccle_copy_number %>%
            filter(copy_number_label != "norm") %>%
            select(dep_map_id, hugo_symbol) %>%
            add_column(is_altered = TRUE)
    },
    {
        ccle_mutations_coding %>%
            select(dep_map_id, hugo_symbol) %>%
            add_column(is_altered = TRUE)
    }
) %>%
    group_by(dep_map_id, hugo_symbol) %>%
    summarise(is_altered = any(is_altered)) %>%
    ungroup() %>%
    unique()


# Tibble of unexpressed genes for each cancer.
unexpressed_gene_tibble <- enframe(confidently_unexpressed_genes) %>%
    dplyr::rename(cancer = "name", hugo_symbol = "value") %>%
    mutate(cancer = str_to_upper(cancer),
           is_unexpressed = TRUE) %>%
    unnest(hugo_symbol)


#### ---- Model data for CRISPR screen ---- ####

check_uniqueness <- function(df) {
    if(any(table(select(df, dep_map_id, hugo_symbol)) > 1)) {
        stop("There are genes in cell lines with multiple effect scores.")
    }
}

# get the density of the values where `essentiality == essentiality_term`
get_density <- function(tib, essentiality_term) {
    from <- 2
    to <- -2
    tib %>%
        filter(essentiality == !!essentiality_term) %>%
        pull(gene_effect) %>%
        density(from = from, to = to)
}

# get the point of overlap of the essential and nonessential distributions
get_essential_gene_effect_overlap <- function(tib, ...) {
    ed <- get_density(tib, "essential")
    ned <- get_density(tib, "nonessential")
    idx <- (ed$y < ned$y) & (dplyr::between(ed$x, -0.9, -0.1))
    poi <- min(ed$x[idx])
    return(poi)
}

cache("model_data_base",
      depends = c("kras_mutation_tib",
                  "essentiality_tib",
                  "cell_lines",
                  "gene_effect",
                  "ccle_expression",
                  "ccle_copy_number",
                  "ccle_mutations_coding"),
{
    info(logger, "Beginning creation of `model_date_base`.")

    # construct base data frame of data for modeling
    model_data_base <- gene_effect %>%
        filter(!dep_map_id %in% c(kras_double_mutants, mapk_muts)) %>%
        filter(dep_map_id %in% !!sample_info_tomerge$dep_map_id) %>%
        left_join(sample_info_tomerge, by = "dep_map_id") %>%
        left_join(kras_mut_tib_tomerge, by = "dep_map_id") %>%
        mutate(
            codon = ifelse(is.na(codon), "WT", codon),
            allele = ifelse(is.na(allele), "WT", allele),
            copy_number_label = ifelse(is.na(copy_number_label), "norm", copy_number_label)
        ) %>%
        mutate(essentiality = case_when(
            hugo_symbol %in% !!essential_genes ~ "essential",
            hugo_symbol %in% !!nonessential_genes ~ "nonessential"
        )) %>%
        filter(!is.na(gene_effect))

    check_uniqueness(model_data_base)

    log_rows(logger, model_data_base, "model_data_base (1)")


#### ---- Filtering alleles and genes from model data ---- ####

    # cutoff 1: 1.5 s.d. from mean of essential genes (per cancer)
    info(logger, "Finding cutoff #1 - 1.5 s.d. from mean of essential genes.")
    cutoff_tib <- model_data_base %>%
        filter(essentiality == "essential") %>%
        group_by(cancer) %>%
        summarise(
            mean_effect = mean(gene_effect),
            sd_effect = sd(gene_effect)
        ) %>%
        mutate(cuttoff_at_2sd = mean_effect + (1.5 * sd_effect)) %>%
        select(cancer, cuttoff_at_2sd)


    # cutoff 2: overlap of dist. of essentials and nonessentials
    info(logger, "Finding cutoff #2 - where dist. of essential and nonessential overlap.")
    cutoff_tib <- model_data_base %>%
        group_by(cancer) %>%
        nest() %>%
        mutate(cutoff_between_non_essential = purrr::map_dbl(data, get_essential_gene_effect_overlap)) %>%
        select(cancer, cutoff_between_non_essential) %>%
        right_join(cutoff_tib, by = "cancer")

    assign("cutoff_tib", cutoff_tib, envir = .GlobalEnv)
    info(logger, "Caching `cutoff_tib`.")
    cache("cutoff_tib")

    info(logger, "Merging the cutoff data frame with the model data.")
    model_data_base <- left_join(model_data_base, cutoff_tib, by = "cancer")

    check_uniqueness(model_data_base)

    log_rows(logger, model_data_base, "model_data_base (2)")


#### ---- Add RNA expression ---- ####

    info(logger, "Adding RNA expression levels for each gene to model data frame.")
    model_data_base <- left_join(
        model_data_base, ccle_expression, by = c("dep_map_id", "hugo_symbol")
    )

    check_uniqueness(model_data_base)


#### ---- Add if target gene is mutated/CNA ---- ####

    info(logger, "Adding column for if gene is altered to model data frame.")
    model_data_base <- left_join(
            model_data_base, gene_is_aletered_tib,
            by = c("dep_map_id", "hugo_symbol")
        ) %>%
        mutate(is_altered = ifelse(is.na(is_altered), FALSE, TRUE))

    check_uniqueness(model_data_base)

    log_rows(logger, model_data_base, "model_data_base (3)")
    info(logger, "Caching `model_data_base`.")


#### ---- Is the gene expressed? ---- ####
    # add a column `is_unexpressed` with T/F

    model_data_base <- model_data_base %>%
        left_join(unexpressed_gene_tibble, by = c("cancer", "hugo_symbol")) %>%
        mutate(is_unexpressed = ifelse(is.na(is_unexpressed), FALSE, TRUE))

    check_uniqueness(model_data_base)
    info(logger, "Caching CRISPR screen model data base.")
    return(model_data_base)
})


cache("model_data", depends = "model_data_base",
{
    # only allele with at least 3 cell lines
    model_data <- model_data_base %>%
        filter(!is_unexpressed) %>%
        group_by(cancer, allele) %>%
        filter(n_distinct(dep_map_id) >= 3) %>%
        ungroup() %>%
        unique()


    # Check that there is no duplicate data.
    if (any(table(model_data$dep_map_id, model_data$hugo_symbol) > 1)) {
        stop("Some data is duplicated in `model_data`.")
        model_data %>%
            group_by(dep_map_id, hugo_symbol) %>%
            filter(n() > 1) %>%
            ungroup()
    }

    info(logger, "Caching `model_data`.")
    return(model_data)
})





#### ---- RNAi Screen Data Prep. ---- ####


cache("rnai_model_data_base",
      depends = c("kras_mutation_tib",
                  "essentiality_tib",
                  "cell_lines",
                  "rnai_effect",
                  "ccle_expression",
                  "ccle_copy_number",
                  "ccle_mutations_coding"),
{
    info(logger, "Beginning creation of `model_date_base`.")

    # construct base data frame of data for modeling
    rnai_model_data_base <- rnai_effect %>%
        filter(!is.na(gene_effect)) %>%
        filter(!dep_map_id %in% c(kras_double_mutants, mapk_muts)) %>%
        filter(dep_map_id %in% !!sample_info_tomerge$dep_map_id) %>%
        left_join(sample_info_tomerge, by = "dep_map_id") %>%
        left_join(kras_mut_tib_tomerge, by = "dep_map_id") %>%
        mutate(
            codon = ifelse(is.na(codon), "WT", codon),
            allele = ifelse(is.na(allele), "WT", allele),
            copy_number_label = ifelse(is.na(copy_number_label), "norm", copy_number_label)
        ) %>%
        mutate(essentiality = case_when(
            hugo_symbol %in% !!essential_genes ~ "essential",
            hugo_symbol %in% !!nonessential_genes ~ "nonessential"
        )) %>%
        filter(cancer != "MM")

    check_uniqueness(rnai_model_data_base)
    log_rows(logger, rnai_model_data_base, "rnai_model_data_base (1)")


#### ---- Filtering alleles and genes from model data ---- ####

    # cutoff 1: 1.5 s.d. from mean of essential genes (per cancer)
    info(logger, "Finding cutoff #1 - 1.5 s.d. from mean of essential genes.")
    rnai_cutoff_tib <- rnai_model_data_base %>%
        filter(essentiality == "essential") %>%
        group_by(cancer) %>%
        summarise(
            mean_effect = mean(gene_effect),
            sd_effect = sd(gene_effect)
        ) %>%
        mutate(cuttoff_at_2sd = mean_effect + (1.5 * sd_effect)) %>%
        select(cancer, cuttoff_at_2sd)


    # cutoff 2: overlap of dist. of essentials and nonessentials
    info(logger, "Finding cutoff #2 - where dist. of essential and nonessential overlap.")
    rnai_cutoff_tib <- rnai_model_data_base %>%
        group_by(cancer) %>%
        nest() %>%
        mutate(cutoff_between_non_essential = purrr::map_dbl(data, get_essential_gene_effect_overlap)) %>%
        select(cancer, cutoff_between_non_essential) %>%
        right_join(rnai_cutoff_tib, by = "cancer")

    assign("rnai_cutoff_tib", rnai_cutoff_tib, envir = .GlobalEnv)
    info(logger, "Caching `rnai_cutoff_tib`.")
    cache("rnai_cutoff_tib")

    info(logger, "Merging the cutoff data frame with the model data.")
    rnai_model_data_base <- left_join(rnai_model_data_base, rnai_cutoff_tib, by = "cancer")

    check_uniqueness(rnai_model_data_base)
    log_rows(logger, rnai_model_data_base, "rnai_model_data_base (2)")


#### ---- Add RNA expression ---- ####

    info(logger, "Adding RNA expression levels for each gene to model data frame.")
    rnai_model_data_base <- left_join(
        rnai_model_data_base, ccle_expression, by = c("dep_map_id", "hugo_symbol")
    )

    check_uniqueness(rnai_model_data_base)


#### ---- Add if target gene is mutated/CNA ---- ####

    info(logger, "Adding column for if gene is altered to model data frame.")
    rnai_model_data_base <- left_join(
            rnai_model_data_base, gene_is_aletered_tib,
            by = c("dep_map_id", "hugo_symbol")
        ) %>%
        mutate(is_altered = ifelse(is.na(is_altered), FALSE, TRUE))

    check_uniqueness(rnai_model_data_base)
    log_rows(logger, rnai_model_data_base, "rnai_model_data_base (3)")
    info(logger, "Caching `rnai_model_data_base`.")


#### ---- Is the gene expressed? ---- ####
    # add a column `is_unexpressed` with T/F

    rnai_model_data_base <- rnai_model_data_base %>%
        left_join(unexpressed_gene_tibble, by = c("cancer", "hugo_symbol")) %>%
        mutate(is_unexpressed = ifelse(is.na(is_unexpressed), FALSE, TRUE))

    check_uniqueness(rnai_model_data_base)
    return(rnai_model_data_base)
})


cache("rnai_model_data", depends = "rnai_model_data_base",
{
    # only allele with at least 3 cell lines
    rnai_model_data <- rnai_model_data_base %>%
        filter(!is_unexpressed) %>%
        filter(!(cancer == "LUAD" & allele == "G13D")) %>%
        group_by(cancer, allele) %>%
        filter(n_distinct(dep_map_id) >= 3) %>%
        ungroup() %>%
        unique()

    check_uniqueness(rnai_model_data)
    info(logger, "Caching `rnai_model_data`.")
    return(rnai_model_data)
})

# rnai_model_data %>%
#     group_by(cancer, allele) %>%
#     summarise(n = n_distinct(dep_map_id)) %>%
#     arrange(cancer, n)



}