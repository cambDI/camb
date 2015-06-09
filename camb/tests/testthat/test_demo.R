library(stringr)
context("Demo Test Function")

test_that("str_length is number of characters", {
  expect_equal(str_length("a"), 1)
  expect_equal(str_length("ab"), 2)
  expect_equal(str_length("abc"), 3)
})

# see http://r-pkgs.had.co.nz/tests.html for more details on testing

