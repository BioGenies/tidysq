# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")

sq_dna <- construct_sq_dna(str_dna, is_clean = TRUE)
sq_ami <- construct_sq_ami(str_ami, is_clean = FALSE)

seqinr_dna <- lapply(str_dna, function(x)
  seqinr::as.SeqFastadna(seqinr::s2c(x)))
seqinr_ami <- lapply(str_ami, function(x)
  seqinr::as.SeqFastaAA(seqinr::s2c(x)))

# IMPORT ----
test_that("correctly imports seqinr::SeqFastadna", {
  expect_identical(import_sq(seqinr_dna, separate = FALSE)[["sq"]],
                   sq_dna)
})

test_that("correctly imports seqinr::SeqFastaAA", {
  expect_identical(import_sq(seqinr_ami, separate = FALSE)[["sq"]],
                   sq_ami)
})

# EXPORT ----
test_that("correctly exports sq object to seqinr::SeqFastadna", {
  expect_identical(export_sq(sq_dna, "seqinr::SeqFastadna"),
                   seqinr_dna)
})

test_that("correctly exports sq object to seqinr::SeqFastaAA", {
  expect_identical(export_sq(sq_ami, "seqinr::SeqFastaAA"),
                   seqinr_ami)
})
