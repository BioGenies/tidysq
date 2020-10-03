# SETUP ----
sq_dna <- construct_sq_dna(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
                           is_clean = TRUE)
sq_rna <- construct_sq_rna(c("UGGCGG", "ACGGUUUCGUU", "UGGAACG", "GGCUCGACAGACUGC"),
                           is_clean = TRUE)

# AMISQ CLNSQ PROTOTYPE ----
test_that("translate() returns clean amino acid sq object", {
  expect_vector(translate(sq_dna),
                ptype = .construct_sq_ptype("ami", is_clean = TRUE),
                size = vec_size(sq_dna))
  expect_vector(translate(sq_rna),
                ptype = .construct_sq_ptype("ami", is_clean = TRUE),
                size = vec_size(sq_rna))
})
