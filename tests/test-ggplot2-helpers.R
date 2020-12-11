
test_that("The sizes of a plot are decided correctly.", {
  expect_equal(unname(decide_size("small")), c(4, 4))
  expect_equal(unname(decide_size("medium")), c(6, 6))
  expect_equal(unname(decide_size("large")), c(8, 8))
  expect_equal(unname(decide_size("tall")), c(4, 8))
  expect_equal(unname(decide_size("wide")), c(8, 4))

  expect_error(decide_size())

  expect_equal(decide_size(width = 5), c(5, 5))
  expect_equal(decide_size(width = 0), c(0, 0))


  expect_equal(decide_size(height = 5), c(5, 5))
  expect_equal(decide_size(height = 0), c(0, 0))
})
