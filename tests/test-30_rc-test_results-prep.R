
test_that("There are no NA values in mutation count columns of `cancer_mut_counts`.", {
  expect_false(any(is.na(cancer_mut_counts$num_samples_per_cancer)))
  expect_false(any(is.na(cancer_mut_counts$num_mut_per_cancer)))
})


test_that("There are no NA values in mutation count columns of `allele_mut_counts`.", {
  expect_false(any(is.na(cancer_mut_counts$num_samples_per_cancer_allele)))
  expect_false(any(is.na(cancer_mut_counts$num_mut_per_cancer_allele)))
})


test_that("There are no NA values in mutation count columns of RC-test results.", {
  expect_false(any(is.na(rc_test_results$num_samples_per_cancer)))
  expect_false(any(is.na(rc_test_results$num_samples_per_cancer_allele)))
  expect_false(any(is.na(rc_test_results$num_mut_per_cancer)))
  expect_false(any(is.na(rc_test_results$num_mut_per_cancer_allele)))
})
