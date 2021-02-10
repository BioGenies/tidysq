# SETUP ----
str_dna_1 <- c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA")
str_dna_2 <- c("ATCTTGAAG", "CATATGCGCTA", "ACGTGTCGA")
str_ami <- "OFUFQBODBUDAFNDBZGFG"
str_unt_1 <- c("!!NV!!XFD!", "P!OQ!!-FI", "SP!!I-F-XXS!")
str_unt_2 <- c("OVNU!!OK!!J", "GOK!MI!N!BB!", "DPOFIN!!")
str_unt_3 <- c("!IF", "", "VF!")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

alph_atp <- c("mA", "mY", "nbA", "nsA")

sq_dna_1 <- sq(str_dna_1, alphabet = "dna_bsc")
sq_dna_2 <- sq(str_dna_2, alphabet = "dna_bsc")
sq_ami <- sq(str_ami, alphabet = "ami_ext")
sq_unt_1 <- sq(str_unt_1, alphabet = "unt", NA_letter = "!")
sq_unt_2 <- sq(str_unt_2, alphabet = "unt", NA_letter = "!")
sq_unt_3 <- sq(str_unt_3, alphabet = "unt", NA_letter = "!")
sq_atp <- sq(str_atp, alphabet = alph_atp)
sq_empty <- sq(character(), alphabet = "rna_bsc")

# PROTOTYPE PRESERVATION ----
test_that("paste() finds common prototype for all arguments (recycling length 1 vectors if necessary)", {
  expect_vector(paste(sq_dna_1, sq_dna_2),
                ptype = vec_ptype_common(sq_dna_1, sq_dna_2),
                size = vec_size_common(sq_dna_1, sq_dna_2))
  expect_vector(paste(sq_unt_1, sq_unt_2, sq_unt_3),
                ptype = vec_ptype_common(sq_unt_1, sq_unt_2, sq_unt_3),
                size = vec_size_common(sq_unt_1, sq_unt_2, sq_unt_3))
  expect_vector(paste(sq_dna_1, sq_unt_3, sq_dna_2),
                ptype = vec_ptype_common(sq_dna_1, sq_unt_3, sq_dna_2),
                size = vec_size_common(sq_dna_1, sq_unt_3, sq_dna_2))
  expect_vector(paste(sq_ami, sq_unt_1),
                ptype = vec_ptype_common(sq_ami, sq_unt_1),
                size = vec_size_common(sq_ami, sq_unt_1))
})

# VALUE COMPUTATION ----
test_that("paste() correctly merges sq objects of the same type and alphabet", {
  expect_identical(
    paste(sq_dna_1, sq_dna_2),
    sq(paste0(str_dna_1, str_dna_2), alphabet = "dna_bsc")
  )
  expect_identical(
    paste(sq_atp, sq_atp),
    sq(paste0(str_atp, str_atp), alphabet = alph_atp)
  )
})

test_that("paste() correctly merges sq objects of different types and alphabets", {
  expect_identical(
    paste(sq_unt_1, sq_unt_2, sq_unt_3),
    vec_cast(
      sq(paste0(str_unt_1, str_unt_2, str_unt_3), alphabet = "unt", NA_letter = "!"),
      vec_ptype_common(sq_unt_1, sq_unt_2, sq_unt_3)
    )
  )
  expect_identical(
    paste(sq_dna_1, sq_unt_3, sq_dna_2),
    vec_cast(
      sq(paste0(str_dna_1, str_unt_3, str_dna_2), alphabet = "unt", NA_letter = "!"),
      vec_ptype_common(sq_dna_1, sq_unt_3, sq_dna_2)
    )
  )
  expect_identical(
    paste(sq_ami, sq_unt_1),
    vec_cast(
      sq(paste0(str_ami, str_unt_1), alphabet = "unt", NA_letter = "!"),
      vec_ptype_common(sq_ami, sq_unt_1)
    )
  )
})

# NO CHANGES IF PASSED 1 ARGUMENT ----
test_that("paste() does nothing to one sq object", {
  expect_identical(paste(sq_unt_3),
                   sq_unt_3)
})
