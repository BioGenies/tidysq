# SETUP ----
na_letter <- .get_na_letter()

char_shortened <- c("A", "C", "G")
char_short <- c("A", "C", "G", "T")
char_medium <- c("A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-")
char_long <- c(LETTERS, letters)

alph_shortened <- sq_alphabet(char_shortened)
alph_short <- sq_alphabet(char_short)
alph_medium <- sq_alphabet(char_medium)
alph_long <- sq_alphabet(char_long)

# CONSTRUCTED ALPHABET VALUE ----
test_that("sq_alphabet object is a subclass of character vector", {
  expect_s3_class(alph_short, "character", exact = FALSE)
  expect_s3_class(alph_medium, "character", exact = FALSE)
  expect_s3_class(alph_long, "character", exact = FALSE)
})

# NOTE: in testthat v.3 expect_equivalent will be deprecated
#  and should be replaced with expect_equal(ignore_attr = TRUE)
test_that("sq_alphabet object coerces to original character vector", {
  expect_equivalent(as.character(alph_short), char_short)
  expect_equivalent(as.character(alph_medium), char_medium)
  expect_equivalent(as.character(alph_long), char_long)
})

# PROTOTYPE ACCEPTANCE ----
test_that("sq_alphabet method accepts character vectors of single characters", {
  expect_vector(alph_short,
                ptype = sq_alphabet_ptype(), size = length(char_short))
  expect_vector(alph_medium,
                ptype = sq_alphabet_ptype(), size = length(char_medium))
  expect_vector(alph_long,
                ptype = sq_alphabet_ptype(), size = length(char_long))
})

# NA_CHARACTER EXTRACTION ----
test_that("sq_alphabet object has na_letter attribute", {
  expect_equal(attr(alph_short, "na_letter"), na_letter)
  expect_equal(attr(alph_medium, "na_letter"), na_letter)
  expect_equal(attr(alph_long, "na_letter"), na_letter)
})

test_that("na_letter method extracts na_letter attribute", {
  expect_equal(na_letter(alph_short), na_letter)
  expect_equal(na_letter(alph_medium), na_letter)
  expect_equal(na_letter(alph_long), na_letter)
})

# INDEXING ----
test_that("indexing works as always for indices different from 2^len - 1 value", {
  expect_equal(alph_short[3], char_short[3])
  expect_equal(alph_medium[2:7], char_medium[2:7])
  expect_equal(alph_long[-6:-19], char_long[-6:-19])
  expect_equal(alph_short[727], char_short[727])
})

test_that("2^len - 1 index extracts na_letter attribute", {
  expect_equal(alph_short[7], na_letter)
  expect_equal(alph_medium[31], na_letter)
  expect_equal(alph_long[63], na_letter)
})

# REMOVING CHARACTERS ----
test_that(".skip_characters remove chosen letters while preserving attributes", {
  expect_equal(
    .skip_characters(alph_short, "T"),
    alph_shortened
  )
  expect_equal(
    .skip_characters(alph_medium, c("U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "-")),
    alph_shortened
  )
})

test_that(".skip_characters doesn't fail while removing letters that aren't present", {
  expect_equal(
    .skip_characters(alph_short, c("T", "X", "Z")),
    alph_shortened
  )
  expect_equal(
    .skip_characters(alph_long, "&"),
    alph_long
  )
})

# EDGE CASES BEHAVIOUR ----
test_that("sq_alphabet method accepts arbitrary strings as characters", {
  expect_vector(sq_alphabet(""),
                ptype = sq_alphabet_ptype(), size = 1)
  expect_vector(sq_alphabet(c("Lorem", "ipsum", "dolor", "sit", "amet")),
                ptype = sq_alphabet_ptype(), size = 5)
})

test_that("sq_alphabet method accepts character vector of length 0", {
  expect_vector(alph_empty <- sq_alphabet(character()),
                ptype = sq_alphabet_ptype(), size = 0)
  expect_s3_class(alph_empty, "character", exact = FALSE)
  expect_equal(na_letter(alph_empty), na_letter)
  expect_equal(alph_empty[1], na_letter)
})
