# Make lollipop plots for the mutations on any gene in any subset of samples.

GRAPHS_DIR <- "20_55_lillipop_plots"
reset_graph_directory(GRAPHS_DIR)

# The cached data frame `cancer_coding_muts_maf` has the proper formatting
# for becoming a `maftools::maf` object.

get_maf_dataframe_for_cancer_allele <- function(cancer, kras_allele) {
    df <- cancer_coding_muts_maf %>%
        filter(cancer %in% !!cancer & ras_allele %in% !!kras_allele)
    return(df)
}

get_maf_dataframe_for_cancer_butallele <- function(cancer, but_kras_allele) {
    df <- cancer_coding_muts_maf %>%
        filter(cancer %in% !!cancer & !(ras_allele %in% !!but_kras_allele))
    return(df)
}

# Make and save a lollipop plot.
lollipop_plot_wrapper <- function(maf_df,
                                  gene,
                                  save_name,
                                  size = c(7, 3)) {
    if (length(size) != 2) stop("`size` must be a vector of length 2.")

    maf <- maftools::read.maf(maf_df, verbose = FALSE)

    svg(plot_path(GRAPHS_DIR, save_name), width = size[[1]], height = size[[2]])
    maftools::lollipopPlot(maf = maf, gene = gene, AACol = "Protein_Change")
    dev.off()
}


# Make and save a lollipop plot.
lollipop2_plot_wrapper <- function(maf_df1,
                                   maf_df2,
                                   gene,
                                   save_name,
                                   top_name = NULL,
                                   bottom_name = NULL,
                                   size = c(9, 3.5)) {
    if (length(size) != 2) stop("`size` must be a vector of length 2.")

    maf1 <- maftools::read.maf(maf_df1, verbose = FALSE)
    maf2 <- maftools::read.maf(maf_df2, verbose = FALSE)

    svg(plot_path(GRAPHS_DIR, save_name), width = size[[1]], height = size[[2]])
    maftools::lollipopPlot2(m1 = maf1, m2 = maf2,
                            gene = gene,
                            AACol1 = "Protein_Change",
                            AACol2 = "Protein_Change",
                            m1_name = top_name,
                            m2_name =  bottom_name)
    dev.off()
}

lollipop_cancer_allele <- function(cancer, kras_allele, gene,
                                   save_template = "lollipop_{cancer}_{str_remove(kras_allele, 'KRAS_')}_{gene}.svg") {
    df <- get_maf_dataframe_for_cancer_allele(cancer, kras_allele)

    save_name <- glue(save_template)
    lollipop_plot_wrapper(df, gene, save_name)
}

lollipop2_cancer_alleles <- function(cancer,
                                    kras_allele1,
                                    kras_allele2,
                                    gene,
                                    save_template = "lollipop_{cancer}_{str_remove(kras_allele1, 'KRAS_')}_{str_remove(kras_allele2, 'KRAS_')}_{gene}.svg") {
    df1 <- get_maf_dataframe_for_cancer_allele(cancer, kras_allele1)
    df2 <- get_maf_dataframe_for_cancer_allele(cancer, kras_allele2)

    save_name <- glue(save_template)
    lollipop2_plot_wrapper(maf_df1 = df1,
                           maf_df2 = df2,
                           gene = gene,
                           save_name = save_name,
                           top_name = str_remove(kras_allele1, "KRAS_"),
                           bottom_name = str_remove(kras_allele2, "KRAS_"))
}

lollipop2_cancer_allele <- function(cancer,
                                    kras_allele,
                                    gene,
                                    save_template = "lollipop_{cancer}_{str_remove(kras_allele, 'KRAS_')}_REST_{gene}.svg") {
    df1 <- get_maf_dataframe_for_cancer_allele(cancer, kras_allele)
    df2 <- get_maf_dataframe_for_cancer_butallele(cancer, kras_allele)

    allele <- str_remove(kras_allele, "KRAS_")
    save_name <- glue(save_template)
    lollipop2_plot_wrapper(maf_df1 = df1,
                           maf_df2 = df2,
                           gene = gene,
                           save_name = save_name,
                           top_name = allele,
                           bottom_name = "rest")
}



tibble::tribble(
    ~cancer, ~kras_allele, ~gene,
    "LUAD", "G13C", "NF1"
) %>%
    mutate(kras_allele = paste0("KRAS_", kras_allele)) %>%
    pwalk(lollipop_cancer_allele)




tibble::tribble(
    ~cancer, ~kras_allele1, ~kras_allele2, ~gene,
    "COAD", "G12D", "G13D", "TP53",
    "COAD", "G12D", "G12V", "APC"
) %>%
    mutate(kras_allele1 = paste0("KRAS_", kras_allele1),
           kras_allele2 = paste0("KRAS_", kras_allele2)) %>%
    pwalk(lollipop2_cancer_alleles)


tibble::tribble(
    ~cancer, ~kras_allele, ~gene,
    "COAD", "A146T", "PORCN",
    "COAD", "G13D", "RICTOR",
    "COAD", "G13D", "RELB",
    "COAD", "G12C", "RAPGEF2",
    "COAD", "G12V", "SMAD4",
    "LUAD", "G13C", "NF1",
    "LUAD", "G12C", "STK11"
) %>%
    mutate(kras_allele = paste0("KRAS_", kras_allele)) %>%
    pwalk(lollipop2_cancer_allele)
