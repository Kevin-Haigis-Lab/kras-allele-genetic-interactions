# The processes for assessing significance of comutation calculations.

# Values for assessing whether or not an interaction is significant.
#   'p_value': p-value thresholds
#   'mutfreq': the mutation frequency of other gene
#   'comutfreq': the frequency of comutation between the allele and other gene
comut_cutoffs <- list(
  p_value = list(mutex = 0.01, comut = 0.01),
  mutfreq = list(mutex = 0.02, comut = 0.01),
  comut_params = list(min_events = 3, comutfreq = 0.10),
  mutex_params = list(t_AM_min = 10)
)


assess_rc_test_significance <- function(rc_df) {
  rc_df %>%
    mutate(
      is_sig =
        t_AM >= comut_cutoffs$mutex_params$t_AM_min &
        num_mut_per_cancer / num_samples_per_cancer > !!comut_cutoffs$mutfreq$mutex &
        p_val < !!comut_cutoffs$p_value$mutex
    )
}


assess_fisher_test_significance <- function(fish_df) {
  fish_df %>%
    mutate(
      is_sig =
        (
          p_value_great < comut_cutoffs$p_value$mutex |
          p_value_less < comut_cutoffs$p_value$comut
        ) &
        n11 >= comut_cutoffs$comut_params$min_events &
        (
          ((n10 + n11) / (n00 + n10 + n01 + n11) > !!comut_cutoffs$mutfreq$comut) |
            (n11 / (n01 + n11) > !!comut_cutoffs$comut_params$comutfreq)
        )
    )
}
