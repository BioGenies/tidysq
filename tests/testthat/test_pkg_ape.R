# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
matrix_dna <- matrix(c("A", "C", "A", "T", "T", "A", "A", "N", "T", "A", "A", "T"), ncol = 3)

sq_dna <- construct_sq_dna(str_dna, is_clean = TRUE)
sq_ami <- construct_sq_ami(str_ami, is_clean = FALSE)
sq_matrix <- construct_sq_dna(apply(matrix_dna, 1, paste, collapse = ""), is_clean = FALSE)

ape_dna <- ape::as.DNAbin(unname(vapply(str_dna, strsplit, list(character()), split = "")))
ape_ami <- ape::as.AAbin(unname(vapply(str_ami, strsplit, list(character()), split = "")))
ape_matrix <- ape::as.DNAbin(matrix_dna)

# IMPORT ----
test_that("correctly imports ape::DNAbin when it's a list", {
  expect_identical(import_sq(ape_dna)[["sq"]],
                   sq_dna)
})

test_that("correctly imports ape::AAbin when it's a list", {
  expect_identical(import_sq(ape_ami)[["sq"]],
                   sq_ami)
})

test_that("correctly imports ape::DNAbin when it's a matrix", {
  expect_identical(import_sq(ape_matrix)[["sq"]],
                   sq_matrix)
})

# EXPORT ----
test_that("correctly exports sq object to ape::DNAbin", {
  expect_identical(export_sq(sq_dna, "ape::DNAbin"),
                   ape_dna)
})

test_that("correctly exports sq object to ape::AAbin", {
  expect_identical(export_sq(sq_ami, "ape::AAbin"),
                   ape_ami)
})
