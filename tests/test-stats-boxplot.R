test_that("p-values are formatted properly.", {
  expect_equal(format_pvalue_label(0.001), "p=1e-03")
  expect_equal(format_pvalue_label(0.00143521), "p=1.44e-03")

  expect_equal(format_pvalue_label(0.001, with_p = FALSE), "1e-03")
  expect_equal(format_pvalue_label(0.00143521, with_p = FALSE), "1.44e-03")
})