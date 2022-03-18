# ERROR FOR NON-SQ OBJECTS ----
test_sq_only <- function(fun, ...) {
  test_that("the first argument must be of sq class", {
    expect_error(fun(1:7, ...))
    expect_error(fun(LETTERS, ...))
    expect_error(fun(list(mean, sum, sd), ...))
  })
}
