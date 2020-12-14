# SETUP ----
str_dna <- c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA")
sq_dna <- sq(str_dna, alphabet = "dna_bsc")

sq_warning <- "Detected letters that do not match specified type!"
bite_warning <- "some sequences are subsetted with index bigger than length - NA introduced"

# WARNING OPTION ----
test_that("on_warning = \"warning\" raises a warning", {
  expect_warning(
    sq(str_dna, "rna_bsc", safe_mode = TRUE, on_warning = "warning"),
    sq_warning
  )
  expect_warning(
    bite(sq_dna, 13:16, on_warning = "warning"),
    bite_warning
  )
})

# ERROR OPTION ----
test_that("on_warning = \"error\" throws an error", {
  expect_error(
    sq(str_dna, "rna_bsc", safe_mode = TRUE, on_warning = "error"),
    sq_warning
  )
  expect_error(
    bite(sq_dna, 13:16, on_warning = "error"),
    bite_warning
  )
})

# MESSAGE OPTION ----
test_that("on_warning = \"message\" produces a message", {
  expect_message(
    sq(str_dna, "rna_bsc", safe_mode = TRUE, on_warning = "message"),
    sq_warning
  )

  skip("warning handling on C++ side is poorly implemented")
  expect_message(
    bite(sq_dna, 13:16, on_warning = "message"),
    bite_warning
  )
})

# SILENT OPTION ----
# TODO: issue #65
test_that("on_warning = \"silent\" suppresses all warnings", {
  expect_silent(sq(str_dna, "rna_bsc", safe_mode = TRUE, on_warning = "silent"))
  expect_silent(bite(sq_dna, 13:16, on_warning = "silent"))
})

# DEFAULT OPTION ----
test_that("default warning handling is to display message as warning", {
  withr::local_options(list(tidysq_on_warning = "warning"))
  expect_warning(
    sq(str_dna, "rna_bsc", safe_mode = TRUE, on_warning = "warning"),
    sq_warning
  )
  expect_warning(
    bite(sq_dna, 13:16, on_warning = "warning"),
    bite_warning
  )
})
