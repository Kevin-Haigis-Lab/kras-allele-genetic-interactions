

test_that("Data is pulled from a ggplot2 object", {
  p <- ggplot(diamonds, aes(x = carat, y = cut)) +
    geom_point()
  expect_equal(pull_original_plot_data(p), diamonds)
})



test_that("Data is pulled from a wrapped ggplot2 object.", {
  p1 <- ggplot(diamonds, aes(x = caret, y = cut)) +
    geom_point()
  p2 <- ggplot(diamonds, aes(x = depth, y = price)) +
    geom_col()
  p <- patchwork::wrap_plots(p1, p2)

  tbl <- pull_wrapped_plot_data(p, number, c("1", "2"))

  expected_tbl <- bind_rows(
    diamonds %>% mutate(number = "1"),
    diamonds %>% mutate(number = "2")
  )

  expect_equal(tbl, expected_tbl)
})



test_that("The correct source data file path is created", {
  expect_null(get_source_data_filename(1))
  expect_equal(
    get_source_data_filename(42, "a"),
    file.path(SOURCE_DATA_BASE_DIR, "figure-01", "panel-a.tsv")
  )
  expect_equal(
    get_source_data_filename(42),
    file.path(SOURCE_DATA_BASE_DIR, "figure-01", "data.tsv")
  )
})



test_that("The directory name is extracted from a file path.", {
  expect_equal(get_directory("a/b/c/d.tsv"), "a/b/c/")
})
