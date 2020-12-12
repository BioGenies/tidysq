# SETUP ----
str_dna_bsc <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_dna_ext <- c("NARTYVTCY", "", "ATKCYGDD", "", "DNAKYTD")
str_rna_bsc <- c("UAUCAGU-A-GU-CA", "CUG-A-CUGAG-CC", "-CUG-AGAGUA-")
str_rna_ext <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_ami_bsc <- c("ACEH", "PASAI", "MALACCA", "SIAK")
str_ami_ext <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
str_unt <- c("VIP01", "VIP02", "VIP04", "MISSING_ONE")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

alph_atp <- c("mA", "mY", "nbA", "nsA")

sq_dna_bsc <- sq(str_dna_bsc, alphabet = "dna_bsc")
sq_dna_ext <- sq(str_dna_ext, alphabet = "dna_ext")
sq_rna_bsc <- sq(str_rna_bsc, alphabet = "rna_bsc")
sq_rna_ext <- sq(str_rna_ext, alphabet = "rna_ext")
sq_ami_bsc <- sq(str_ami_bsc, alphabet = "ami_bsc")
sq_ami_ext <- sq(str_ami_ext, alphabet = "ami_ext")
sq_unt <- sq(str_unt, alphabet = "unt")
sq_atp <- sq(str_atp, alphabet = alph_atp)

# SQ TYPE IS STRONGER THAN CHARACTER VECTOR ----
# Testing common prototypes as below was impossible with pure c() due to how
# R dispatches methods and how we cannot change dispatch from character class.
test_that("character vector and sq object have common ptype that is of the same sq type as sq", {
  expect_identical(vec_ptype2(str_dna_bsc, sq_dna_bsc),
                   vec_ptype(sq_dna_bsc))
  expect_identical(vec_ptype2(str_unt, sq_dna_ext),
                   vec_ptype(sq_dna_ext))
  expect_identical(vec_ptype2(str_rna_ext, sq_rna_bsc),
                   vec_ptype(sq_rna_bsc))
  expect_identical(vec_ptype2(str_rna_bsc, sq_rna_ext),
                   vec_ptype(sq_rna_ext))
  expect_identical(vec_ptype2(str_ami_bsc, sq_ami_bsc),
                   vec_ptype(sq_ami_bsc))
  expect_identical(vec_ptype2(str_atp, sq_ami_ext),
                   vec_ptype(sq_ami_ext))
  expect_identical(vec_ptype2(str_dna_ext, sq_unt),
                   vec_ptype(sq_unt))
  expect_identical(vec_ptype2(str_ami_ext, sq_atp),
                   vec_ptype(sq_atp))
})
