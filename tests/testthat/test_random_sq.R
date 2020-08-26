# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("random_sq() returns an object of sq class", {
  expect_vector(random_sq(4, 7, "dna", TRUE),
                ptype = .construct_sq_ptype("dna", is_clean = TRUE),
                size = 4)
  expect_vector(random_sq(1, 26, "ami", FALSE),
                ptype = .construct_sq_ptype("ami", is_clean = FALSE),
                size = 1)
})

test_that("each sequence of random_sq() is of the passed length", {
  for (sq in random_sq(3, 27, "rna", TRUE)) {
    expect_equal(attr(sq, "original_length"), 27)
  }
})

# LENGTH SAFETY ----
test_that("using sd argument of random_sq() doesn't generate negative-length sequences", {
  # Using local mock to ensure that rnorm generates some negative values
  # Will have to change it if random_sq() changes length generating function
  # This might suggest that we should extract this as a separate function
  local_mock(rnorm = function(x, mean, sd) {
    seq(from = mean - 2*sd, to = mean + 2*sd, length.out = x)
  })
  for (sq in random_sq(25, 3, "dna", TRUE, sd = 50)) {
    expect_gte(attr(sq, "original_length"), 0)
  }
})

# EDGE CASES ----
test_that("random_sq() can generate 0 sequences", {
  expect_vector(random_sq(0, 13, "rna", FALSE),
                ptype = .construct_sq_ptype("rna", is_clean = FALSE),
                size = 0)
})

test_that("random_sq() can generate 0-length sequences", {
  empty_sq <- random_sq(3, 0, "ami", TRUE)
  expect_vector(empty_sq,
                ptype = .construct_sq_ptype("ami", is_clean = TRUE),
                size = 3)
  for (sq in empty_sq) {
    expect_length(sq, 0)
  }
})
