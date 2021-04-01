# SETUP ----
sq_dna <- sq(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
             alphabet = "dna_bsc")
sq_rna <- sq(c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-"),
             alphabet = "rna_ext")
sq_ami <- sq(c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR"),
             alphabet = "ami_ext")
sq_unt <- sq(c("vip01", "vip02", "vip04", "missing_one"),
             alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""), alphabet = c("mA", "mY", "nbA", "nsA"))

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("sq_type() returns valid sq type", {
  expect_sq_type(sq_type(sq_dna), unt.ok = TRUE, atp.ok = TRUE)
  expect_sq_type(sq_type(sq_rna), unt.ok = TRUE, atp.ok = TRUE)
  expect_sq_type(sq_type(sq_ami), unt.ok = TRUE, atp.ok = TRUE)
  expect_sq_type(sq_type(sq_unt), unt.ok = TRUE, atp.ok = TRUE)
  expect_sq_type(sq_type(sq_atp), unt.ok = TRUE, atp.ok = TRUE)
})

# ERROR FOR NON-SQ OBJECTS ----
test_that("sq_type() throws an error whenever passed object of class other that sq", {
  expect_error(sq_type(1:7))
  expect_error(sq_type(LETTERS))
  expect_error(sq_type(list(mean, sum, sd)))
})

test_that("`sq_type<-`() throws an error whenever passed object of class other that sq", {
  expect_error(sq_type(1:7) <- "ami_bsc")
  expect_error(sq_type(LETTERS) <- "curiosity killed the cat")
  expect_error(sq_type(list(mean, sum, sd)) <- "wdym")
})

# CORRECT VALUE ON ACCESS ----
test_that("sq_type() returns correct sq type", {
  expect_equal(sq_type(sq_dna),
               "dna_bsc")
  expect_equal(sq_type(sq_rna),
               "rna_ext")
  expect_equal(sq_type(sq_ami),
               "ami_ext")
  expect_equal(sq_type(sq_unt),
               "unt")
  expect_equal(sq_type(sq_atp),
               "atp")
})

# TYPIFY ON ASSIGN ----
test_that("`sq_type<-`() is an alias for typify", {
  expect_identical(
    {
      sq_type(sq_dna) <- "ami_bsc"
      sq_dna
    },
    typify(sq_dna, "ami_bsc")
  )
  expect_identical(
    {
      sq_type(sq_rna) <- "ami_ext"
      sq_rna
    },
    typify(sq_rna, "ami_ext")
  )
})
