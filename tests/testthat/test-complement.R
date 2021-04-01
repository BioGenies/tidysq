# SETUP ----
sq_dna <- sq(c("ATCTTGAAG", "CATATGCGCTA", "ACGTGTCGA", ""),
             alphabet = "dna_bsc")
sq_dna_compl <- sq(c("TAGAACTTC", "GTATACGCGAT", "TGCACAGCT", ""),
                   alphabet = "dna_bsc")
sq_dna_2 <- sq(c("KCYSRRCACNB", "BAYRNYWAK", "NBVKAWRYGG"),
               alphabet = "dna_ext")
sq_dna_2_compl <- sq(c("MGRSYYGTGNV", "VTRYNRWTM", "NVBMTWYRCC"),
                     alphabet = "dna_ext")
sq_rna <- sq(c("UAGUAACCGUAAGCG", "UAGUCC--UA-G"),
             alphabet = "rna_bsc")
sq_rna_compl <- sq(c("AUCAUUGGCAUUCGC", "AUCAGG--AU-C"),
                   alphabet = "rna_bsc")

# PROTOTYPE PRESERVATION ----
test_that("complement() preserves all attributes of original vector", {
  expect_vector(complement(sq_dna),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(complement(sq_rna),
                ptype = vec_ptype(sq_rna),
                size = vec_size(sq_rna))
  expect_vector(complement(sq_dna_2),
                ptype = vec_ptype(sq_dna_2),
                size = vec_size(sq_dna_2))
})

# ERROR FOR NON-DNA/RNA OBJECTS ----
test_that("complement() throws an error whenever passed object of class other that sq_dna/sq_rna", {
  expect_error(complement(19:8))
  expect_error(complement(list(mean, sum, sd)))
  expect_error(complement(LETTERS))
  expect_error(complement(sq(character(), "ami_bsc")))
  expect_error(complement(sq(c("accmsce", "auprcacc"), alphabet = c("auprc", "acc", "msce"))))
})

# VALUE COMPUTATION ----
# NOTE: used as.character() because hypothetically one value might have
#  multiple equivalent representations
test_that("complement() returns correct complement value for complement-only characters", {
  expect_equal(
    as.character(complement(sq_dna)),
    as.character(sq_dna_compl)
  )
})

test_that("complement() returns correct complement value for complement-less characters", {
  expect_equal(
    as.character(complement(sq_rna)),
    as.character(sq_rna_compl)
  )
})

test_that("complement() returns correct complement value for ambiguous characters", {
  expect_equal(
    as.character(complement(sq_dna_2)),
    as.character(sq_dna_2_compl)
  )
})

# CANCELLING UPON DOUBLE USAGE ----
test_that("double use of complement() returns original value", {
  expect_identical(complement(complement(sq_dna)), sq_dna)
  expect_identical(complement(complement(sq_rna)), sq_rna)
  expect_identical(complement(complement(sq_dna_2)), sq_dna_2)
})
