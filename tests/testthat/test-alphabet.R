# SETUP ----
NA_letter <- getOption("tidysq_NA_letter")

char_shortened <- c("A", "C", "G")
char_short <- c("A", "C", "G", "T")
char_medium <- c("A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-")
char_long <- c(LETTERS, letters)

alph_shortened <- sq_alphabet(char_shortened, "atp")
alph_short <- sq_alphabet(char_short, "atp")
alph_medium <- sq_alphabet(char_medium, "rna_ext")
alph_long <- sq_alphabet(char_long, "atp")

sq_dna <- sq(c("ACTGTC", "CGCGTTA"), alphabet = "dna_bsc")
sq_ami <- sq(c("APOPNIQEV", "CSVMIBF"), alphabet = "ami_ext")
sq_unt <- sq(c("GO%NC@E(123)RO", "NFI%(#)VT;"), alphabet = "unt")

# CONSTRUCTED ALPHABET VALUE ----
test_that("sq_alphabet object is a subclass of character vector", {
  expect_s3_class(alph_short, "character", exact = FALSE)
  expect_s3_class(alph_medium, "character", exact = FALSE)
  expect_s3_class(alph_long, "character", exact = FALSE)
})

test_that("sq_alphabet object coerces to original character vector", {
  expect_equal(as.character(alph_short), char_short, ignore_attr = TRUE)
  expect_equal(as.character(alph_medium), char_medium, ignore_attr = TRUE)
  expect_equal(as.character(alph_long), char_long, ignore_attr = TRUE)
})

# ALPHABET EXTRACTION ----
test_that("get_sq_alphabet() extracts an object of sq_alphabet class", {
  expect_vector(alphabet(sq_dna),
                ptype = sq_alphabet_ptype("dna_bsc"))
  expect_vector(alphabet(sq_ami),
                ptype = sq_alphabet_ptype("ami_ext"))
  expect_vector(alphabet(sq_unt),
                ptype = sq_alphabet_ptype("unt"))
})

test_that("get_sq_alphabet() extracts \"alphabet\" attribute", {
  expect_identical(alphabet(sq_dna),
                   attr(sq_dna, "alphabet"))
  expect_identical(alphabet(sq_ami),
                   attr(sq_ami, "alphabet"))
  expect_identical(alphabet(sq_unt),
                   attr(sq_unt, "alphabet"))
})

# PROTOTYPE ACCEPTANCE ----
test_that("sq_alphabet method accepts character vectors of single characters", {
  expect_vector(alph_short,
                ptype = sq_alphabet_ptype("atp"),
                size = length(char_short))
  expect_vector(alph_medium,
                ptype = sq_alphabet_ptype("rna_ext"),
                size = length(char_medium))
  expect_vector(alph_long,
                ptype = sq_alphabet_ptype("atp"),
                size = length(char_long))
})


# INDEXING ----
test_that("indexing works as always for indices different from 2^len - 1 value", {
  expect_equal(alph_short[3], char_short[3])
  expect_equal(alph_medium[2:7], char_medium[2:7])
  expect_equal(alph_long[-6:-19], char_long[-6:-19])
  expect_equal(alph_short[727], char_short[727])
})

test_that("2^len - 1 index extracts default NA_letter attribute when not specified", {
  expect_equal(alph_short[7], NA_letter)
  expect_equal(alph_medium[31], NA_letter)
  expect_equal(alph_long[63], NA_letter)
})

test_that("2^len - 1 index extracts specified NA_letter", {
  expect_equal(alph_short[7, NA_letter = "?"], "?")
  expect_equal(alph_medium[31, NA_letter = NA_letter], NA_letter)
  expect_equal(alph_long[63, NA_letter = "NA"], "NA")
})

# EDGE CASES BEHAVIOUR ----
test_that("sq_alphabet method accepts arbitrary strings as characters", {
  expect_vector(sq_alphabet("", "unt"),
                ptype = sq_alphabet_ptype("unt"),
                size = 1)
  expect_vector(sq_alphabet(c("Lorem", "ipsum", "dolor", "sit", "amet"), "atp"),
                ptype = sq_alphabet_ptype("atp"),
                size = 5)
})

test_that("sq_alphabet method accepts character vector of length 0", {
  expect_vector(alph_empty <- sq_alphabet(character(), "unt"),
                ptype = sq_alphabet_ptype("unt"),
                size = 0)
  expect_s3_class(alph_empty, "character", exact = FALSE)
  expect_equal(alph_empty[0], NA_letter)
})
