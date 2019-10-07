
# Prints all of the output file names for Snakemake
# it isn't elegant, but is sufficient

library(dplyr)

# cancer data
cancer_data <- cancer_coding_muts_df %>%
    filter(!is_hypermutant)

make_output_filename <- function(cancer, ras_allele, ...) {
    rasallele <- stringr::str_remove_all(ras_allele, "_")
    fnames <- c()
    for (mut_test in c("exclusivity", "comutation")) {
        fbasename <- paste0(cancer, "_", rasallele, "_", mut_test, "_results.rds")
        fnames <- c(fnames, file.path("output", mut_test, fbasename))
    }
    return(fnames)
}

fnames <- cancer_data %>%
    filter(cancer != "SKCM") %>%
    group_by(cancer, ras_allele) %>%
    summarise(num_samples = n_distinct(case_id)) %>%
    ungroup() %>%
    filter(num_samples >= 10 & ras_allele != "WT") %>%
    filter(!stringr::str_detect(ras_allele, "nonfssub")) %>%
    purrr::pmap(make_output_filename) %>%
    unlist()

for (fname in fnames) {
    cat("\"", fname, "\"", ",\n", sep = "")
}

# "output/exclusivity/COAD_KRASA146T_exclusivity_results.rds",
# "output/comutation/COAD_KRASA146T_comutation_results.rds",
# "output/exclusivity/COAD_KRASG12A_exclusivity_results.rds",
# "output/comutation/COAD_KRASG12A_comutation_results.rds",
# "output/exclusivity/COAD_KRASG12C_exclusivity_results.rds",
# "output/comutation/COAD_KRASG12C_comutation_results.rds",
# "output/exclusivity/COAD_KRASG12D_exclusivity_results.rds",
# "output/comutation/COAD_KRASG12D_comutation_results.rds",
# "output/exclusivity/COAD_KRASG12S_exclusivity_results.rds",
# "output/comutation/COAD_KRASG12S_comutation_results.rds",
# "output/exclusivity/COAD_KRASG12V_exclusivity_results.rds",
# "output/comutation/COAD_KRASG12V_comutation_results.rds",
# "output/exclusivity/COAD_KRASG13D_exclusivity_results.rds",
# "output/comutation/COAD_KRASG13D_comutation_results.rds",
# "output/exclusivity/LUAD_KRASG12A_exclusivity_results.rds",
# "output/comutation/LUAD_KRASG12A_comutation_results.rds",
# "output/exclusivity/LUAD_KRASG12C_exclusivity_results.rds",
# "output/comutation/LUAD_KRASG12C_comutation_results.rds",
# "output/exclusivity/LUAD_KRASG12D_exclusivity_results.rds",
# "output/comutation/LUAD_KRASG12D_comutation_results.rds",
# "output/exclusivity/LUAD_KRASG12V_exclusivity_results.rds",
# "output/comutation/LUAD_KRASG12V_comutation_results.rds",
# "output/exclusivity/LUAD_KRASG13C_exclusivity_results.rds",
# "output/comutation/LUAD_KRASG13C_comutation_results.rds",
# "output/exclusivity/MM_KRASG12A_exclusivity_results.rds",
# "output/comutation/MM_KRASG12A_comutation_results.rds",
# "output/exclusivity/MM_KRASG12D_exclusivity_results.rds",
# "output/comutation/MM_KRASG12D_comutation_results.rds",
# "output/exclusivity/MM_KRASG12R_exclusivity_results.rds",
# "output/comutation/MM_KRASG12R_comutation_results.rds",
# "output/exclusivity/MM_KRASG12V_exclusivity_results.rds",
# "output/comutation/MM_KRASG12V_comutation_results.rds",
# "output/exclusivity/MM_KRASG13D_exclusivity_results.rds",
# "output/comutation/MM_KRASG13D_comutation_results.rds",
# "output/exclusivity/MM_KRASQ61H_exclusivity_results.rds",
# "output/comutation/MM_KRASQ61H_comutation_results.rds",
# "output/exclusivity/MM_KRASQ61L_exclusivity_results.rds",
# "output/comutation/MM_KRASQ61L_comutation_results.rds",
# "output/exclusivity/MM_KRASQ61R_exclusivity_results.rds",
# "output/comutation/MM_KRASQ61R_comutation_results.rds",
# "output/exclusivity/PAAD_KRASG12C_exclusivity_results.rds",
# "output/comutation/PAAD_KRASG12C_comutation_results.rds",
# "output/exclusivity/PAAD_KRASG12D_exclusivity_results.rds",
# "output/comutation/PAAD_KRASG12D_comutation_results.rds",
# "output/exclusivity/PAAD_KRASG12R_exclusivity_results.rds",
# "output/comutation/PAAD_KRASG12R_comutation_results.rds",
# "output/exclusivity/PAAD_KRASG12V_exclusivity_results.rds",
# "output/comutation/PAAD_KRASG12V_comutation_results.rds",
# "output/exclusivity/PAAD_KRASQ61H_exclusivity_results.rds",
# "output/comutation/PAAD_KRASQ61H_comutation_results.rds",
# "output/exclusivity/PAAD_KRASQ61R_exclusivity_results.rds",
# "output/comutation/PAAD_KRASQ61R_comutation_results.rds"
