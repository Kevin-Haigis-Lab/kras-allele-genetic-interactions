
files_to_source <- c(
    "10_11_syn-let_heatmaps-boxplots.R",
    "10_15_linear-modeling-syn-let_ppi-subnetworks.R",
    "10_37_gsea-depmap-analysis.R",
    "10_45_overall-dependencies-by-allele.R",
    "10_55_paad_depmap_jun-cdkn2a-G12V.R",
    "10_57_coad_depmap_wdr26.R",
    "20_40_highlivel-genetic-interactions.R",
    "20_41_disagreeing-interactions_logOR-barplot.R",
    "20_43_apriori-lists-genetic-interactions.R",
    "20_45_fxnal-enrich-genetic-interactions.R",
    "20_47_enriched-functions_signaling-pathways.R",
    "20_48_enriched-functions_compare-functions_heatmaps.R",
    "20_50_rainfall-plots.R",
    "40_15_comparing-COAD-allele-subnetworks.R",
    "40_17_comparing-PAAD-allele-subnetworks.R",
    "40_20_comut-dependency-genes-ppi-connectivity.R",
    "60_10_MM-specific-oncogenes.R",
    "70_10_survival-analysis.R",
    "70_15_comutation-survival-analysis.R",
    "70_30_survival-analysis-hits-followup.R",
    "90_03_mutation-burden-distribution.R",
    "90_05_kras-allele-distribution.R",
    "90_10_NF1-depletion-COAD.R"
)


for (f in files_to_source) {
    cat(glue("running: '{f}'"), "\n")
    source(file.path("src", f))
}


build_comutation_figure(c(22:27, 29:35, 38))