test_that("helper_functions: check_list", {
  list_to_check <- c("a"=1, "b"=2, "c"=3, "d"=4, "e"=5)
  expect_equal(check_list(list_to_check, c("a", "e")), TRUE)
  expect_equal(check_list(list_to_check, "c"), TRUE)
  expect_equal(check_list(list_to_check, "v"), FALSE)
  expect_equal(check_list(list_to_check, c("a", "v")), FALSE)
  expect_equal(check_list(list_to_check, ""), FALSE)
})
