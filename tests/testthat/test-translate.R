# SETUP ----
sq_dna <- sq(c("TACTGGGCATGA", "CAGGTCGGA", "TAGTAGTCC", "ACGGTG"),
             alphabet = "dna_bsc")
sq_rna <- sq(c("UGGCGG", "ACGGUUUCGUU", "UGGAACG", "GGCUCGACAGACUGC", ""),
             alphabet = "rna_bsc")

str_dna_translated <- c("YWA*", "QVG", "**S", "TV")
str_rna_translated <- c("WR", "TVS", "WN", "GSTDC", "")

# SQ_AMI_BSC PROTOTYPE ----
test_that("translate() returns clean amino acid sq object", {
  expect_vector(translate(sq_dna),
                ptype = sq_ptype(get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_dna))
  expect_vector(translate(sq_rna),
                ptype = sq_ptype(get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_rna))
})

# VALUE COMPUTATION ----
test_that("translate() returns correct value", {
  expect_equivalent(as.character(translate(sq_dna)),
                    str_dna_translated)
  expect_equivalent(as.character(translate(sq_rna)),
                    str_rna_translated)
})
