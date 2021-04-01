# SETUP ----
sq_1 <- sq(c("QWERTYUIOP", "ASDF-GHJKL", "ZXCV-BNM"), alphabet = "unt")
sq_2 <- sq(c("", "CAGTGT", "CGGCTATXT"), alphabet = LETTERS)
sq_3 <- sq(c("AreYouOK", "WhoAreYou", "YesWeCan"),
           alphabet = c("Are", "You", "Who", "Yes", "We", "Can", "O", "K"))
sq_4 <- sq(character(), "dna_bsc")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("find_invalid_letters() returns a list of correct length", {
  expect_vector(find_invalid_letters(sq_1, "dna_ext"),
                ptype = list(),
                size = vec_size(sq_1))
})

test_that("each element of find_invalid_letters() is a character vector", {
  for (inv_letters in find_invalid_letters(sq_1, "dna_ext")) {
    expect_vector(inv_letters,
                  ptype = character())
  }
})

# ERROR FOR NON-SQ OBJECTS ----
test_that("find_invalid_letters() throws an error whenever passed object of class other that sq", {
  expect_error(find_invalid_letters(1:7, "dna_bsc"))
  expect_error(find_invalid_letters(LETTERS, "ami_ext"))
  expect_error(find_invalid_letters(list(mean, sum, sd), "rna_ext"))
})

# VALUE COMPUTATION ----
test_that("find_invalid_letters() correctly computes value", {
  expect_equal(
    find_invalid_letters(sq_1, "dna_ext"),
    list(c("E", "I", "O", "P", "Q", "U"), c("F", "J", "L"), c("X", "Z"))
  )
})

test_that("find_invalid_letters() consider NA part of alphabet", {
  expect_equal(
    find_invalid_letters(sq_2, "dna_bsc"),
    list(character(), character(), "X")
  )
})

test_that("find_invalid_letters() works for alphabet with multiple-character letters", {
  expect_equal(
    find_invalid_letters(sq_3, "ami_ext"),
    list(c("Are", "You"), c("Are", "You", "Who"), c("Yes", "We", "Can"))
  )
})

test_that("find_invalid_letters() returns value without additional attributes", {
  expect_equal(
    find_invalid_letters(sq_1, "dna_ext"),
    list(c("E", "I", "O", "P", "Q", "U"), c("F", "J", "L"), c("X", "Z"))
  )
  expect_equal(
    find_invalid_letters(sq_2, "dna_bsc"),
    list(character(), character(), "X")
  )
  expect_equal(
    find_invalid_letters(sq_3, "ami_ext"),
    list(c("Are", "You"), c("Are", "You", "Who"), c("Yes", "We", "Can"))
  )
})

# EDGE CASES ----
test_that("passing empty sq to find_invalid_letters() returns empty list()", {
  expect_equal(
    find_invalid_letters(sq_4, "rna_bsc"),
    list()
  )
})
