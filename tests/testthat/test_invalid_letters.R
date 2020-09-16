# SETUP ----
sq_1 <- construct_sq(c("QWERTYUIOP", "ASDF-GHJKL", "ZXCV-BNM"))
sq_2 <- construct_sq(c("", "CAGTgT", "CGGCTAtXT"),
                     non_standard = LETTERS)
sq_3 <- construct_sq(c("AreYouOK", "WhoAreYou", "YesWeCan"),
                     non_standard = c("Are", "You", "Who", "Yes", "We", "Can", "O", "K"))
sq_4 <- construct_sq(character(), "dna")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("get_invalid_letters() returns a list of correct length", {
  expect_vector(get_invalid_letters(sq_1, "dna"),
                ptype = list(),
                size = vec_size(sq_1))
})

test_that("each element of get_invalid_letters() is a character vector", {
  for (inv_letters in get_invalid_letters(sq_1, "dna")) {
    expect_vector(inv_letters,
                  ptype = character())
  }
})

# VALUE COMPUTATION ----
test_that("get_invalid_letters() correctly computes value", {
  expect_equivalent(
    get_invalid_letters(sq_1, "dna"),
    list(c("Q", "E", "U", "I", "O", "P"), c("F", "J", "L"), c("Z", "X"))
  )
})

test_that("get_invalid_letters() consider NA part of alphabet", {
  expect_equivalent(
    get_invalid_letters(sq_2, "dna"),
    list(character(), character(), "X")
  )
})

test_that("get_invalid_letters() works for alphabet with multiple-character letters", {
  expect_equivalent(
    get_invalid_letters(sq_3, "ami"),
    list(c("Are", "You"), c("Who", "Are", "You"), c("Yes", "We", "Can"))
  )
})

test_that("get_invalid_letters() returns value without additional attributes", {
  expect_equal(
    get_invalid_letters(sq_1, "dna"),
    list(c("Q", "E", "U", "I", "O", "P"), c("F", "J", "L"), c("Z", "X"))
  )
  expect_equal(
    get_invalid_letters(sq_2, "dna"),
    list(character(), character(), "X")
  )
  expect_equal(
    get_invalid_letters(sq_3, "ami"),
    list(c("Are", "You"), c("Who", "Are", "You"), c("Yes", "We", "Can"))
  )
})

# EDGE CASES ----
test_that("passing empty sq to get_invalid_letters() returns empty list()", {
  expect_equal(
    get_invalid_letters(sq_4, "rna"),
    list()
  )
})
