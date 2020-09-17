# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_dna_rev <- c("GTACGGGTCAT", "AGGCTGGAC", "GCCTGATGAT", "", "TGGCA")
str_ami <- c("PADINI", "HOUDINI", "UNBELIEVABLE")
str_ami_rev <- c("INIDAP", "INIDUOH", "ELBAVEILEBNU")

sq_dna <- construct_sq_dna(str_dna, is_clean = TRUE)
sq_ami <- construct_sq_ami(str_ami, is_clean = FALSE)

# PROTOTYPE PRESERVATION ----
test_that("reverse() preserves all attributes of original vector", {
  expect_vector(reverse(sq_dna),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(reverse(sq_ami),
                ptype = vec_ptype(sq_ami),
                size = vec_size(sq_ami))
})

# VALUE COMPUTATION ----
test_that("reverse() returns correct value", {
  expect_equivalent(as.character(reverse(sq_dna)),
                    str_dna_rev)
  expect_equivalent(as.character(reverse(sq_ami)),
                    str_ami_rev)
})

# CANCELLING UPON DOUBLE USAGE ----
test_that("double use of reverse() returns original value", {
  expect_identical(reverse(reverse(sq_dna)),
                   sq_dna)
  expect_identical(reverse(reverse(sq_ami)),
                   sq_ami)
})
