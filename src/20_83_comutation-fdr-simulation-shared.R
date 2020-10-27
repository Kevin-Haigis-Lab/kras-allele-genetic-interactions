
temp_sim_dir <- file.path("temp", "comut-fdr")
temp_input_dir <- file.path(temp_sim_dir, "input")
temp_output_dir <- file.path(temp_sim_dir, "output")

simulation_input_file_name <- function(idx) {
  file.path(temp_input_dir, glue("simulation_{idx}.tsv"))
}

simulation_output_file_name <- function(idx) {
  file.path(temp_output_dir, glue("simulation_{idx}.qs"))
}


# Ensure there is always a 2x2 contingency table.
my_table <- function(a, b) {
  a <- factor(a, levels = c("FALSE", "TRUE"))
  b <- factor(b, levels = c("FALSE", "TRUE"))
  table(a, b)
}


# Would the simulation pass the project's standard for comutation?
project_significance_check <- function(mut_table, p_value) {

  ## THE FOLLOWING PARAMETERS WERE COPIED FROM "20_35_rc-fisher-comparison.R"
  # p-value thresholds for mutual exclusivity and comutation
  p_val_cut_mutex <- 0.01
  p_val_cut_comut <- 0.01

  # Thresholds for the mutation frequency of other gene for mut. ex. and comut.
  mutfreq_mutex <- 0.02
  mutfreq_comut <- 0.01

  # Threshold for the frequency of comutation between the allele and other gene.
  comutfreq_comut <- 0.10
  ############################################################################

  if (p_value >= p_val_cut_comut) {
    return(FALSE)
  }
  if (mut_table[2, 2] == 0) {
    return(FALSE)
  }

  num_samples_per_cancer <- sum(mut_table)
  mut_freq_1 <- sum(mut_table[2, ]) / num_samples_per_cancer
  mut_freq_2 <- sum(mut_table[, 2]) / num_samples_per_cancer

  if (mut_freq_1 == 0 | mut_freq_2 == 0) {
    return(FALSE)
  }

  comut_freq_1 <- mut_table[2, 2] / sum(mut_table[2, ])
  comut_freq_2 <- mut_table[2, 2] / sum(mut_table[, 2])

  check_1 <- (mut_freq_1 > mutfreq_comut | comut_freq_1 > comutfreq_comut)
  check_2 <- (mut_freq_2 > mutfreq_comut | comut_freq_2 > comutfreq_comut)

  return(check_1 & check_2)
}
