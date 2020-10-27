#!/bin/env Rscript

## Run all of the analyses for the comutation paper.

library(ProjectTemplate)
load.project()

run_src_files <- function(files) {
  for (f in files) {
    src(file.path("src", f))
  }
}


#### ---- Run all analysese scripts ---- ####
message("Starting to run analyses scripts.")

src_files <- c(
  "05_10_expression-filter-plots.R",
  "05_20_possible-hotspot-mutations.R",
  "10_05_describe-cell-lines.R",
  "10_10_linear-modeling-syn-let.R",
  "10_11_syn-let_heatmaps-boxplots.R",
  "10_13_linear-modeling-syn-let_fxnal-enrichment.R",
  "10_15_linear-modeling-syn-let_ppi-subnetworks.R",
  "10_17_linear-modeling-syn-let_specific-protein-types.R",
  "10_18_alluvial-plots-comparing-dependency.R",
  "10_30_gsea-depmap-preparation.R"
)
run_src_files(src_files)

system("./src/10_34_submit-gsea-depmap.sh")
readline(prompt = "Hit enter when all GSEA runs have completed.")

src_files <- c(
  "10_37_gsea-depmap-analysis.R",
  "10_45_overall-dependencies-by-allele.R",
  "10_50_single-sample-rank-sum.R",
  "10_55_paad_depmap_jun-cdkn2a-G12V.R",
  "10_56_luad_depmap_G12C-cisplatin-sensitivity.R",
  "10_57_coad_depmap_wdr26.R",
  "20_02_prepare-input-tables.R"
)
run_src_files(src_files)

system("./src/20_21_run-rc-test-snakemake.sh")
readline(prompt = "Hit enter when the Snakemake workflow has completed.")

src_files <- c(
  "20_29_rc-test_results-prep.R",
  "20_30_rc-test-results-analysis.R",
  "20_34_rc-fisher-assessment-processes.R",
  "20_35_rc-fisher-comparison.R",
  "20_40_highlivel-genetic-interactions.R",
  "20_41_disagreeing-interactions_logOR-barplot.R",
  "20_43_apriori-lists-genetic-interactions.R",
  "20_44_select-enriched-functions.R",
  "20_45_fxnal-enrich-genetic-interactions.R",
  "20_46_enriched-functions_bar-plots.R",
  "20_47_enriched-functions_signaling-pathways.R",
  "20_48_enriched-functions_compare-functions_heatmaps.R",
  "20_49_rainfall-plots-subroutines.R",
  "20_50_rainfall-plots.R",
  "20_55_lollipop_plots.R",
  "20_60_nonallelespec_prepare-input-tables.R"
)
run_src_files(src_files)


system("./src/20_63_nonallelespec_run-rc-test-snakemake.sh")
readline(prompt = "Hit enter when the Snakemake workflow has completed.")

src_files <- c(
  "20_65_nonallele-specific-comutation-results-preparation.R",
  "20_66_nonallele-specific-increased-comutation.R",
  "20_67_non-allele-specific-comutation-comparison.R",
  "20_70_luad-g12c-stk11.R",
  "20_75_allele-to-allele-comutation-analysis.R",
  "30_10_logisitic-regression-syn-let.R",
  "30_20_discriminant-analysis-syn-let.R",
  "40_05_describe-combined-ppi.R",
  "40_10_overlap-synlet-comutation.R",
  "40_12_overlap-synlet-comutation-hits.R",
  "40_15_comparing-COAD-allele-subnetworks.R",
  "40_16_comparing-LUAD-allele-subnetworks.R",
  "40_17_comparing-PAAD-allele-subnetworks.R",
  "40_20_comut-dependency-genes-ppi-connectivity.R",
  "40_25_luad_co-enriched-functions.R",
  "40_50_comparing-alleles-across-cancer.R",
  "40_60_synlet-explained-by-comuts.R",
  "40_62_genes-for-synlet-comut.R",
  "40_63_explore-synlet-comut.R",
  "40_64_synlet-comut_describe-muts.R",
  "40_65_overview-plot-synlet-comut.R",
  "50_09_example_observed-predicted-kras-alleles.R",
  "50_10_observed-predicted-kras-alleles.R",
  "50_11_observed-predicted-kras-alleles_v2.R",
  "50_12_observed-predicted-kras-alleles_v3.R",
  "50_13_observed-predicted-kras-alleles_v3_per_tsb.R",
  "50_15_correlate-mutsig-allele-probs.R",
  "50_20_mutsignatures-distributions.R",
  "50_30_mutsignatures_prob-causing-allele.R",
  "50_35_mutational-signature-allele-associations.R",
  "50_60_example-mutational-signature-spectra.R",
  "60_10_MM-specific-oncogenes.R",
  "60_20_mutyh-signature18-associations.R",
  "70_10_survival-analysis.R",
  "70_15_comutation-survival-analysis.R",
  "70_30_survival-analysis-hits-followup.R",
  "90_03_mutation-burden-distribution.R",
  "90_05_kras-allele-distribution.R",
  "90_10_NF1-depletion-COAD.R",
  "90_13_comutation-network_withingroup-ppi-closeness.R",
  "90_20_dataset-description-table.R",
  "90_34_random-description-plots.R"
)
run_src_files(src_files)

message("Analyses complete")

#### ---- Prepare paper assets ---- ####
message("Compiling paper assets.")

compile_supp_data()
build_comutation_figures()
copy_final_figures()

message("Done")
