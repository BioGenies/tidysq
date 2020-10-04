# SETUP ----
sq_dna <- construct_sq_dna(c("TACTGGGCATGA", "CAGGTCGGA", "TAGTAGTCC", "ACGGTG"),
                           is_clean = TRUE)
sq_rna <- construct_sq_rna(c("UGGCGG", "ACGGUUUCGUU", "UGGAACG", "GGCUCGACAGACUGC", ""),
                           is_clean = TRUE)

str_dna_translated <- c("YWA*", "QVG", "**S", "TV")
str_rna_translated <- c("WR", "TVS", "WN", "GSTDC", "")

# AMISQ CLNSQ PROTOTYPE ----
test_that("translate() returns clean amino acid sq object", {
  expect_vector(translate(sq_dna),
                ptype = .construct_sq_ptype("ami", is_clean = TRUE),
                size = vec_size(sq_dna))
  expect_vector(translate(sq_rna),
                ptype = .construct_sq_ptype("ami", is_clean = TRUE),
                size = vec_size(sq_rna))
})

# VALUE COMPUTATION ----
test_that("translate() returns correct value", {
  expect_equivalent(as.character(translate(sq_dna)),
                    str_dna_translated)
  expect_equivalent(as.character(translate(sq_rna)),
                    str_rna_translated)
})
