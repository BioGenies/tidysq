# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_ami <- c("REGENERATED", "TECHNICAL", "FEAT")
matrix_dna <- matrix(c("A", "C", "A", "T", "T", "A", "A", "G", "T", "A", "A", "T"), ncol = 3)

sq_dna <- sq(str_dna, alphabet = "dna_bsc")
sq_ami <- sq(str_ami, alphabet = "ami_bsc")
sq_matrix <- sq(apply(matrix_dna, 1, paste, collapse = ""), alphabet = "dna_bsc")

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
