# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_1_dna <- "TCYYCAHGGCHA"
str_rna <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_1_rna <- "UUACGACGGCGCAU"
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
str_1_ami <- "PASVKNFYD"
str_unt <- c("vip01", "vip02", "vip04", "missing_one")
str_1_unt <- "2high4U"

sq_dna <- construct_sq_dna(str_dna, is_clean = TRUE)
sq_1_dna <- construct_sq_dna(str_1_dna, is_clean = FALSE)
sq_rna <- construct_sq_rna(str_rna, is_clean = FALSE)
sq_1_rna <- construct_sq_rna(str_1_rna, is_clean = TRUE)
sq_ami <- construct_sq_ami(str_ami, is_clean = FALSE)
sq_1_ami <- construct_sq_ami(str_1_ami, is_clean = TRUE)
sq_unt <- construct_sq(str_unt, type = "unt")
sq_1_unt <- construct_sq(str_1_unt, type = "unt")

biostr_dna <- Biostrings::DNAStringSet(str_dna)
biostr_1_dna <- Biostrings::DNAString(str_1_dna)
biostr_rna <- Biostrings::RNAStringSet(str_rna)
biostr_1_rna <- Biostrings::RNAString(str_1_rna)
biostr_ami <- Biostrings::AAStringSet(str_ami)
biostr_1_ami <- Biostrings::AAString(str_1_ami)
biostr_unt <- Biostrings::BStringSet(str_unt)
biostr_1_unt <- Biostrings::BString(str_1_unt)

# IMPORT ----
test_that("correctly imports Biostrings::DNAStringSet", {
  expect_identical(import_sq(biostr_dna)[["sq"]],
                   sq_dna)
})
test_that("correctly imports Biostrings::DNAString", {
  expect_identical(import_sq(biostr_1_dna)[["sq"]],
                   sq_1_dna)
})

test_that("correctly imports Biostrings::RNAStringSet", {
  expect_identical(import_sq(biostr_rna)[["sq"]],
                   sq_rna)
})
test_that("correctly imports Biostrings::RNAString", {
  expect_identical(import_sq(biostr_1_rna)[["sq"]],
                   sq_1_rna)
})

test_that("correctly imports Biostrings::AAStringSet", {
  expect_identical(import_sq(biostr_ami)[["sq"]],
                   sq_ami)
})
test_that("correctly imports Biostrings::AAString", {
  expect_identical(import_sq(biostr_1_ami)[["sq"]],
                   sq_1_ami)
})

test_that("correctly imports Biostrings::BStringSet", {
  expect_identical(import_sq(biostr_unt)[["sq"]],
                   sq_unt)
})
test_that("correctly imports Biostrings::BString", {
  expect_identical(import_sq(biostr_1_unt)[["sq"]],
                   sq_1_unt)
})

# EXPORT ----
test_that("correctly exports sq object to Biostrings::DNAStringSet", {
  expect_identical(export_sq(sq_dna, "Biostrings::DNAStringSet"),
                   biostr_dna)
})
test_that("correctly exports sq object to Biostrings::DNAString", {
  expect_identical(export_sq(sq_1_dna, "Biostrings::DNAString"),
                   biostr_1_dna)
})

test_that("correctly exports sq object to Biostrings::RNAStringSet", {
  expect_identical(export_sq(sq_rna, "Biostrings::RNAStringSet"),
                   biostr_rna)
})
test_that("correctly exports sq object to Biostrings::RNAString", {
  expect_identical(export_sq(sq_1_rna, "Biostrings::RNAString"),
                   biostr_1_rna)
})

test_that("correctly exports sq object to Biostrings::AAStringSet", {
  expect_identical(export_sq(sq_ami, "Biostrings::AAStringSet"),
                   biostr_ami)
})
test_that("correctly exports sq object to Biostrings::AAString", {
  expect_identical(export_sq(sq_1_ami, "Biostrings::AAString"),
                   biostr_1_ami)
})
