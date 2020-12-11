# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_ami <- c("REGENERATED", "TECHNICAL", "FEAT")

sq_dna <- sq(str_dna, alphabet = "dna_bsc")
sq_ami <- sq(str_ami, alphabet = "ami_bsc")

seqinr_dna <- lapply(str_dna, function(x)
  seqinr::as.SeqFastadna(seqinr::s2c(x)))
seqinr_ami <- lapply(str_ami, function(x)
  seqinr::as.SeqFastaAA(seqinr::s2c(x)))

# Names below are just chosen for extravagance, not political reasons
names_1 <- c("balalaika", "perestroyka", "sojuz", "oncePutin", "LenOut")
names_2 <- c("Unionists", "Confederates", "Mexicans")

seqinr_dna_n <- mapply(
  function(x, name) seqinr::as.SeqFastadna(seqinr::s2c(x), name = name),
  str_dna, names_1, SIMPLIFY = FALSE, USE.NAMES = FALSE
)
seqinr_ami_n <- mapply(
  function(x, name) seqinr::as.SeqFastaAA(seqinr::s2c(x), name = name),
  str_ami, names_2, SIMPLIFY = FALSE, USE.NAMES = FALSE
)

# IMPORT ----
test_that("correctly imports seqinr::SeqFastadna", {
  expect_identical(import_sq(seqinr_dna, separate = FALSE)[["sq"]],
                   sq_dna)
  expect_identical(import_sq(seqinr_dna_n, separate = FALSE)[["sq"]],
                   sq_dna)
  expect_identical(import_sq(seqinr_dna_n, separate = FALSE)[["name"]],
                   names_1)
})

test_that("correctly imports seqinr::SeqFastaAA", {
  expect_identical(import_sq(seqinr_ami, separate = FALSE)[["sq"]],
                   sq_ami)
  expect_identical(import_sq(seqinr_ami_n, separate = FALSE)[["sq"]],
                   sq_ami)
  expect_identical(import_sq(seqinr_ami_n, separate = FALSE)[["name"]],
                   names_2)
})

# EXPORT ----
test_that("correctly exports sq object to seqinr::SeqFastadna", {
  expect_identical(export_sq(sq_dna, "seqinr::SeqFastadna"),
                   seqinr_dna)
  expect_identical(export_sq(sq_dna, "seqinr::SeqFastadna", name = names_1),
                   seqinr_dna_n)
})

test_that("correctly exports sq object to seqinr::SeqFastaAA", {
  expect_identical(export_sq(sq_ami, "seqinr::SeqFastaAA"),
                   seqinr_ami)
  expect_identical(export_sq(sq_ami, "seqinr::SeqFastaAA", name = names_2),
                   seqinr_ami_n)
})
