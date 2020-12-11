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

# CASTING TO CHARACTER ----
test_that("vec_cast() to character works like as.character()", {
  expect_identical(vec_cast(sq_dna_bsc, character()),
                   as.character(sq_dna_bsc))
  expect_identical(vec_cast(sq_dna_ext, character()),
                   as.character(sq_dna_ext))
  expect_identical(vec_cast(sq_rna_bsc, character()),
                   as.character(sq_rna_bsc))
  expect_identical(vec_cast(sq_rna_ext, character()),
                   as.character(sq_rna_ext))
  expect_identical(vec_cast(sq_ami_bsc, character()),
                   as.character(sq_ami_bsc))
  expect_identical(vec_cast(sq_ami_ext, character()),
                   as.character(sq_ami_ext))
  expect_identical(vec_cast(sq_unt, character()),
                   as.character(sq_unt))
  expect_identical(vec_cast(sq_atp, character()),
                   as.character(sq_atp))
})

# CASTING TO STANDARD SQ TYPE ----
test_that("unt sq objects can be casted with vec_cast() to standard sq types", {
  expect_identical(
    vec_cast(sq(str_dna_bsc, alphabet = "unt"), vec_ptype(sq_dna_bsc)),
    sq_dna_bsc
  )
  expect_identical(
    vec_cast(sq(str_dna_ext, alphabet = "unt"), vec_ptype(sq_dna_ext)),
    sq_dna_ext
  )
  expect_identical(
    vec_cast(sq(str_rna_bsc, alphabet = "unt"), vec_ptype(sq_rna_bsc)),
    sq_rna_bsc
  )
  expect_identical(
    vec_cast(sq(str_rna_ext, alphabet = "unt"), vec_ptype(sq_rna_ext)),
    sq_rna_ext
  )
  expect_identical(
    vec_cast(sq(str_ami_bsc, alphabet = "unt"), vec_ptype(sq_ami_bsc)),
    sq_ami_bsc
  )
  expect_identical(
    vec_cast(sq(str_ami_ext, alphabet = "unt"), vec_ptype(sq_ami_ext)),
    sq_ami_ext
  )
})
