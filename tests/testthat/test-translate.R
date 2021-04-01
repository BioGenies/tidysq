# SETUP ----
sq_dna <- sq(c("TACTGGGCATGA", "CAGGTCGGA", "TAGTAGTCC", "ACGGTG"),
             alphabet = "dna_bsc")
sq_dna_2 <- sq(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
               alphabet = "dna_bsc")
sq_rna <- sq(c("UGGCGG", "ACGGUUUCGUU", "UGGAACG", "GGCUCGACAGACUGC", ""),
             alphabet = "rna_bsc")

str_dna_translated <- c("YWA*", "QVG", "**S", "TV")
str_dna_translated_4 <- c("YWAW", "QVG", "**S", "TV")
str_dna_2_translated_31 <- c("YWA", "QVG", "EES", "", "T")
str_dna_2_translated_31_as_stop <- c("YWA", "QVG", "**S", "", "T")
str_dna_translated_12 <- c("YWA*", "QVG", "**S", "TV")
str_rna_translated <- c("WR", "TVS", "WN", "GSTDC", "")

# SQ_AMI_BSC PROTOTYPE ----
test_that("translate() returns clean amino acid sq object", {
  expect_vector(translate(sq_dna),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_dna))
  expect_vector(translate(sq_rna),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_rna))
  expect_vector(translate(sq_dna, table = 12),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_dna))
})

# ERROR FOR NON-DNA/RNA OBJECTS ----
test_that("translate() throws an error whenever passed object of class other that sq_dna/sq_rna", {
  expect_error(translate(19:8))
  expect_error(translate(list(mean, sum, sd)))
  expect_error(translate(LETTERS))
  expect_error(translate(sq(character(), "ami_bsc")))
  expect_error(translate(sq(c("accmsce", "auprcacc"), alphabet = c("auprc", "acc", "msce"))))
})

# VALUE COMPUTATION ----
test_that("translate() returns correct value", {
  expect_equal(as.character(translate(sq_dna)),
               str_dna_translated)
  expect_equal(as.character(translate(sq_rna)),
               str_rna_translated)
})

test_that("translate() correctly handles tables other than default 1", {
  expect_equal(as.character(translate(sq_dna, 4)),
               str_dna_translated_4)
  expect_equal(as.character(translate(sq_dna, 12)),
               str_dna_translated_12)
  expect_equal(as.character(translate(sq_rna, 24)),
               str_rna_translated)
})
