
test_that("The correct word length is retrieved.", {
  l <- expect_message(get_string_length("hi"), "number of chars")
  expect_equal(l, 2)

  l <- expect_message(get_string_length("WT"), "known string widths")
  expect_true(dplyr::near(l, 1))

  expect_equal(ncol(string_length_dictionary), 2)
  expect_equal(colnames(string_length_dictionary), c("word", "known_length"))
})