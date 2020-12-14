
# test `get_hugo_from_depmap_ids()`
test_that("Hugo symbols are extracted.", {
  test_genes <- c(
    "ABCG8 (64241)",
    "ACCSL (390110)",
    "ACTL7A (10881)",
    "ACTL7B (10880)",
    "ACTL9 (284382)",
    "ADAD1 (132612)",
    "ADAM18 (8749)",
    "ADAM2 (2515)",
    "ADAM20 (8748)",
    "ADAM30 (11085)"
  )

  expected_output <- c(
    "ABCG8",
    "ACCSL",
    "ACTL7A",
    "ACTL7B",
    "ACTL9",
    "ADAD1",
    "ADAM18",
    "ADAM2",
    "ADAM20",
    "ADAM30"
  )

  expect_length(get_hugo_from_depmap_ids(test_genes), length(test_genes))
  expect_equal(get_hugo_from_depmap_ids(test_genes), expected_output)
  expect_length(get_hugo_from_depmap_ids(NULL), 0)
})

# test `get_entrez_from_depmap_ids()`
test_that("Entrez IDs are extracted.", {
  test_genes <- c(
    "ABCG8 (64241)",
    "ACCSL (390110)",
    "ACTL7A (10881)",
    "ACTL7B (10880)",
    "ACTL9 (284382)",
    "ADAD1 (132612)",
    "ADAM18 (8749)",
    "ADAM2 (2515)",
    "ADAM20 (8748)",
    "ADAM30 (11085)"
  )

  expected_output <- c(
    "64241",
    "390110",
    "10881",
    "10880",
    "284382",
    "132612",
    "8749",
    "2515",
    "8748",
    "11085"
  )

  expect_length(get_entrez_from_depmap_ids(test_genes), length(test_genes))
  expect_equal(get_entrez_from_depmap_ids(test_genes), expected_output)
  expect_equal(
    get_entrez_from_depmap_ids(test_genes, convert_to_num = TRUE),
    as.numeric(expected_output)
  )
  expect_length(get_entrez_from_depmap_ids(NULL), 0)
})


test_that("Some small helper functions work.", {

  # `replace_na_zero()`
  expect_equal(replace_na_zero(NA), 0)
  expect_equal(replace_na_zero(rep(NA, 3)), rep(0, 3))
  expect_equal(replace_na_zero(c(NA, 0, 2, NA, NA, -5, 4)), c(0, 0, 2, 0, 0, -5, 4))


  # `replace_numeric_NAs()`
  df <- tibble(x = c(NA, 1:3, NA, NA, -2))
  nonzero_df <- replace_numeric_NAs(df)
  expect_equal(nonzero_df$x, replace_na_zero(df$x))

  df$y <- c("A", "B", "C", NA, "E", "F", "G")
  nonzero_df <- replace_numeric_NAs(df)
  expect_equal(nonzero_df$x, replace_na_zero(df$x))
  expect_equal(nonzero_df$y, df$y)

  df$z <- df$x
  nonzero_df <- replace_numeric_NAs(df)
  expect_equal(nonzero_df$x, replace_na_zero(df$x))
  expect_equal(nonzero_df$z, replace_na_zero(df$z))
  expect_equal(nonzero_df$y, df$y)


  # `str_rep()`
  expect_equal(str_rep("A", 4), "AAAA")
  expect_equal(str_rep("AB", 4), "ABABABAB")
  expect_equal(str_rep("", 4), "")
  expect_equal(str_rep(" ", 4), "    ")
  expect_equal(str_rep("A", 1), "A")
  expect_equal(str_rep("A", 0), "")


  # `str_round()`
  expect_equal(str_round(1.111, 1), "1.1")
  expect_equal(str_round(1.1, 3), "1.100")
  expect_equal(str_round(1.111, 4), "1.1110")
  expect_equal(str_round(1.1, 1, 3), "1.100")
  expect_equal(str_round(1.123, 1, 3), "1.100")
})

test_that("Citation functions work properly.", {
  expect_true(is.character(list_all_packages()))
  expect_true(length(list_all_packages()) > 100)  # I have a lot of packages...
})
