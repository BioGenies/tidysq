# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("random_sq() returns an object of sq class", {
  expect_vector(random_sq(4, 7, "dna_bsc"),
                ptype = sq_ptype(CPP_get_standard_alphabet("dna_bsc"), "dna_bsc"),
                size = 4)
  expect_vector(random_sq(1, 26, "ami_ext"),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_ext"), "ami_ext"),
                size = 1)
  expect_vector(random_sq(5, 8, c("mA", "nY")),
                ptype = sq_ptype(c("mA", "nY"), "atp"),
                size = 5)
})

test_that("each sequence of random_sq() is of the passed length", {
  for (sq in random_sq(3, 27, "rna_bsc")) {
    expect_equal(attr(sq, "original_length"), 27)
  }
})

# ARGUMENT PREREQUISITES ----
test_that("random_sq() doesn't accept \"unt\" as type", {
  expect_error(random_sq(5, 11, "unt"),
               "method 'random_sq' cannot take 'unt' as alphabet type")
})

# LENGTH SAFETY ----
test_that("using sd argument of random_sq() doesn't generate negative-length sequences", {
  # Using local mock to ensure that rnorm generates some negative values
  # Will have to change it if random_sq() changes length generating function
  # This might suggest that we should extract this as a separate function
  local_mock(rnorm = function(x, mean, sd) {
    seq(from = mean - 2*sd, to = mean + 2*sd, length.out = x)
  })
  for (sq in random_sq(25, 3, "dna_bsc", sd = 50)) {
    expect_gte(attr(sq, "original_length"), 0)
  }
})

# SEED SAFETY ---
test_that("generating random sequences with the same seed gives the same sequences", {
  set.seed(6125)
  sq_1 <- random_sq(10, 100, "ami_bsc")
  sq_2 <- random_sq(5, 20, "dna_ext", sd = 5)
  set.seed(6125)
  expect_equal(random_sq(10, 100, "ami_bsc"), sq_1)
  expect_equal(random_sq(5, 20, "dna_ext", sd = 5), sq_2)
})

# EDGE CASES ----
test_that("random_sq() can generate 0 sequences", {
  expect_vector(random_sq(0, 13, "rna_ext"),
                ptype = sq_ptype(CPP_get_standard_alphabet("rna_ext"), "rna_ext"),
                size = 0)
})

test_that("random_sq() can generate 0-length sequences", {
  empty_sq <- random_sq(3, 0, "ami_bsc")
  expect_vector(empty_sq,
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = 3)
  for (sq in empty_sq) {
    expect_length(sq, 0)
  }
})
