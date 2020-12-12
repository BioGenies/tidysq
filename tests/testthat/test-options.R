test_that("get_tidysq_options() return a named list", {
  expect_list(get_tidysq_options(),
              names = "unique")
})
