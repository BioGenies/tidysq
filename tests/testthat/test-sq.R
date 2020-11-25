# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_dna_ic <- c("TaCTggGcAtg", "cAggTCgGA", "tAGTAgtCCG", "", "acgGT")
str_rna <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
str_unt <- c("vip01", "vip02", "vip04", "missing_one")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

dna_bsc_alph <- c("A", "C", "G", "T", "-")
rna_ext_alph <- c("A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", 
                  "H", "V", "N", "-")

ami_ext_alph <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
                  "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", 
                  "Z", "-", "*")
atp_alph <- c("mA", "mY", "nbA", "nsA")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("sq() returns object of correct prototype", {
  expect_vector(sq(str_dna, alphabet = "dna_bsc"),
                ptype = sq_ptype(dna_bsc_alph, "dna_bsc"),
                size = vec_size(str_dna))
  expect_vector(sq(str_rna, alphabet = "rna_ext"),
                ptype = sq_ptype(rna_ext_alph, "rna_ext"),
                size = vec_size(str_rna))
  expect_vector(sq(str_ami, alphabet = "ami_ext"),
                ptype = sq_ptype(ami_ext_alph, "ami_ext"),
                size = vec_size(str_ami))
})

test_that("sq() returns object of correct class for unt and atp options", {
  expect_s3_class(sq(str_atp, alphabet = atp_alph),
                  class = "sq_atp",
                  exact = FALSE)
  expect_s3_class(sq(str_unt, alphabet = "unt"),
                  class = "sq_unt",
                  exact = FALSE)
})

test_that("sq() returns object of same size as passed character vector", {
  expect_equal(vec_size(sq(str_atp, alphabet = atp_alph)),
               vec_size(str_atp))
  expect_equal(vec_size(sq(str_unt, alphabet = "unt")),
               vec_size(str_unt))
})

test_that("sq() returns object with alphabet attribute that contains existing letters for unt and atp options", {
  expect_setequal(
    alphabet(sq(str_atp, alphabet = atp_alph)),
    atp_alph
  )
  expect_setequal(
    alphabet(sq(str_unt, alphabet = "unt")),
    obtain_alphabet(str_unt)
  )
})

# NA WHEN ACTUAL ALPHABET MISMATCHES ----
test_that("letters not in alphabet are loaded as NA's ", {
  expect_equivalent(
    as.character(sq(str_ami, "rna_bsc", NA_letter = "!"), NA_letter = "!"),
    c("!U!!A!!!!!", "U!!!!UC!U!!!!!", "!!A!")
  )
  expect_equivalent(
    as.character(sq(str_rna, "ami_bsc", NA_letter = "!"), NA_letter = "!"), 
    c("", "K!S-!VW-AWWWG", "YGHHH-", "-CRASH", "MND-K!!!V-MY-")
  )
})

# ALPHABET UNT WHEN SAFE MODE ----
test_that("type set as untyped when in safe mode and alphabet mismatches", {
  expect_warning(
    sq(str_ami, "rna_bsc", safe_mode = TRUE),
    "Detected letters that do not match specified type!"
  )
  suppressWarnings({
    expect_equivalent(
      as.character(sq(str_ami, "rna_bsc", NA_letter = "!", safe_mode = TRUE)),
      str_ami
    )
    expect_equivalent(
      as.character(sq(str_rna, "ami_bsc", NA_letter = "!", safe_mode = TRUE)), 
      str_rna
    )
  })
})

# IGNORE CASE ----
test_that("ignore_case parameter works correctly", {
  expect_equal(
    sq(str_dna, ignore_case = TRUE),
    sq(str_dna, ignore_case = FALSE)
  )
  expect_equal(
    sq(str_dna, "dna_bsc", ignore_case = TRUE),
    sq(str_dna, "dna_bsc", ignore_case = FALSE)
  )
})


# REVERSING TO CHARACTER ----
test_that("applying as.character() returns original character vector", {
  expect_equivalent(
    as.character(sq(str_dna, alphabet = "dna_bsc")),
    str_dna
  )
  expect_equivalent(
    as.character(sq(str_rna, alphabet = "rna_ext")),
    str_rna
  )
  expect_equivalent(
    as.character(sq(str_ami, alphabet = "ami_ext")),
    str_ami
  )
  expect_equivalent(
    as.character(sq(str_atp, alphabet = atp_alph)),
    str_atp
  )
  expect_equivalent(
    as.character(sq(str_unt, alphabet = "unt")),
    str_unt
  )
})

# TYPE GUESSING ----
test_that("sq() correctly guesses sq type", {
  expect_identical(sq(str_dna, alphabet = "dna_bsc"),
                   sq(str_dna))
  expect_identical(sq(str_rna, alphabet = "rna_ext"),
                   sq(str_rna))
  expect_identical(sq(str_ami, alphabet = "ami_ext"),
                   sq(str_ami))
})