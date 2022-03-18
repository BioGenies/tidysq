# ERROR FOR NON-SQ OBJECTS ----
test_sq_only <- function(fun, ..., .data.frame_ok = FALSE) {
  test_message <- if (.data.frame_ok)
    "the first argument must be of sq or data.frame class" else
      "the first argument must be of sq class"
  
  test_that(test_message, {
    expect_s3_class(rlang::catch_cnd(fun(1:7, ...)),
                    "error_no_method")
    expect_s3_class(rlang::catch_cnd(fun(LETTERS, ...)),
                    "error_no_method")
    expect_s3_class(rlang::catch_cnd(fun(list(mean, sum, sd), ...)),
                    "error_no_method")
  })
}
