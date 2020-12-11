# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_ami <- c("REGENERATED", "TECHNICAL", "FEAT")
matrix_dna <- matrix(c("A", "C", "A", "T", "T", "A", "A", "G", "T", "A", "A", "T"), ncol = 3)
matrix_ami <- matrix(c("H", "M", "M", "A", "D", "F", "Y", "G", "R", "K", "K", "P"), ncol = 4)

sq_dna <- sq(str_dna, alphabet = "dna_bsc")
sq_ami <- sq(str_ami, alphabet = "ami_bsc")
sq_dna_mat <- sq(apply(matrix_dna, 1, paste, collapse = ""), alphabet = "dna_bsc")
sq_ami_mat <- sq(apply(matrix_ami, 1, paste, collapse = ""), alphabet = "ami_bsc")

ape_dna <- ape::as.DNAbin(unname(vapply(str_dna, strsplit, list(character()), split = "")))
ape_ami <- ape::as.AAbin(unname(vapply(str_ami, strsplit, list(character()), split = "")))
ape_dna_mat <- ape::as.DNAbin(matrix_dna)
ape_ami_mat <- ape::as.AAbin(matrix_ami)

# Names below are just chosen for extravagance, not political reasons
names_dna <- c("balalaika", "perestroyka", "sojuz", "oncePutin", "LenOut")
names_ami <- c("Unionists", "Confederates", "Mexicans")

ape_dna_n <- ape::as.DNAbin(setNames(
  vapply(str_dna, strsplit, list(character()), split = ""), names_dna))
ape_ami_n <- ape::as.AAbin(setNames(
  vapply(str_ami, strsplit, list(character()), split = ""), names_ami))

# IMPORT ----
test_that("correctly imports ape::DNAbin when it's a list", {
  expect_identical(import_sq(ape_dna)[["sq"]],
                   sq_dna)
  expect_identical(import_sq(ape_dna_n)[["sq"]],
                   sq_dna)
  expect_identical(import_sq(ape_dna_n)[["name"]],
                   names_dna)
})

test_that("correctly imports ape::AAbin when it's a list", {
  expect_identical(import_sq(ape_ami)[["sq"]],
                   sq_ami)
  expect_identical(import_sq(ape_ami_n)[["sq"]],
                   sq_ami)
  expect_identical(import_sq(ape_ami_n)[["name"]],
                   names_ami)
})

test_that("correctly imports ape::DNAbin when it's a matrix", {
  expect_identical(import_sq(ape_dna_mat)[["sq"]],
                   sq_dna_mat)
})

test_that("correctly imports ape::AAbin when it's a matrix", {
  expect_identical(import_sq(ape_ami_mat)[["sq"]],
                   sq_ami_mat)
})

test_that("correctly imports ape::DNAbin when it's a character vector", {
  expect_identical(import_sq(ape_dna[[1]])[["sq"]],
                   sq_dna[1])
})

test_that("correctly imports ape::AAbin when it's a character vector", {
  expect_identical(import_sq(ape_ami[[1]])[["sq"]],
                   sq_ami[1])
})

# EXPORT ----
test_that("correctly exports sq object to ape::DNAbin", {
  expect_identical(export_sq(sq_dna, "ape::DNAbin"),
                   ape_dna)
  expect_identical(export_sq(sq_dna, "ape::DNAbin", name = names_dna),
                   ape_dna_n)
})

test_that("correctly exports sq object to ape::AAbin", {
  expect_identical(export_sq(sq_ami, "ape::AAbin"),
                   ape_ami)
  expect_identical(export_sq(sq_ami, "ape::AAbin", name = names_ami),
                   ape_ami_n)
})
