# SETUP ----
sq_ami <- sq(c("AGNTYIKFGGAYTIB", "MATEGILIAADGYTWIL", "MIPADHICAANGIENAGIK"),
             alphabet = "ami_ext")
sq_dna <- sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", "CAGACCANNNATAG"),
             alphabet = "dna_ext")
sq_rna <- sq(c("GGCUGCGGGACUGAGGC", "UUCAUGCGGCUAGGGCU", "UAGGCGAGCAGAGUAG"),
             alphabet = "rna_bsc")
sq_unt <- sq(c("GO%NC@E(123)RO", "NFI%(#)VT;"), alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
             alphabet = c("mA", "mY", "nbA", "nsA"))

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("%has% returns a logical vector", {
  expect_vector(sq_ami %has% "A",
                ptype = logical(),
                size = vec_size(sq_ami))
  expect_vector(sq_dna %has% "CTG",
                ptype = logical(),
                size = vec_size(sq_dna))
  expect_vector(sq_rna %has% c("GGG", "UGA", "CUGC"),
                ptype = logical(),
                size = vec_size(sq_rna))
  expect_vector(sq_unt %has% c("NFI", ";"),
                ptype = logical(),
                size = vec_size(sq_unt))
  expect_vector(sq_atp %has% "mYmY",
                ptype = logical(),
                size = vec_size(sq_atp))
})

# ERROR FOR NON-SQ OBJECTS ----
test_that("%has% throws an error whenever passed object of class other that sq", {
  expect_error(7:2 %has% "ATC")
  expect_error(LETTERS %has% "H")
  expect_error(list(mean, sum, sd) %has% c("ALFA", "BUICK", "ROMEO"))
})

# VALUE COMPUTATION FOR SINGLE MOTIFS ----
test_that("%has% works correctly for basic letters", {
  expect_equal(sq_ami %has% "A",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_ami %has% "TY",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_dna %has% "C",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_dna %has% "TAA",
               c(TRUE, TRUE, FALSE))
  expect_equal(sq_rna %has% "C",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_rna %has% "GGGA",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_unt %has% "(123)",
               c(TRUE, FALSE))
  expect_equal(sq_atp %has% "mYmY",
               c(TRUE, FALSE, FALSE))
})

test_that("%has% correctly interprets ambiguous letters in a motif", {
  expect_equal(sq_ami %has% "J", 
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_dna %has% "CANA",
               c(FALSE, FALSE, TRUE))
  expect_equal(sq_rna %has% "GBD",
               c(TRUE, TRUE, TRUE))
})

# WORKING ^ AND $ SIGNS ----
test_that("^ at the beginning matches only from the beginning of a motif", {
  expect_equal(sq_ami %has% "^IM",
               c(FALSE, FALSE, FALSE))
  expect_equal(sq_dna %has% "^C",
               c(TRUE, FALSE, TRUE))
  expect_equal(sq_rna %has% "^U",
               c(FALSE, TRUE, TRUE))
  expect_equal(sq_unt %has% "^NFI",
               c(FALSE, TRUE))
  expect_equal(sq_atp %has% "^mA",
               c(TRUE, FALSE, FALSE))
})

test_that("$ at the end matches only to the end of a motif", {
  expect_equal(sq_ami %has% "IB$",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_dna %has% "AT$",
               c(TRUE, TRUE, FALSE))
  expect_equal(sq_rna %has% "A$",
               c(FALSE, FALSE, FALSE))
  expect_equal(sq_unt %has% "(#)$",
               c(FALSE, FALSE))
  expect_equal(sq_atp %has% "mA$",
               c(FALSE, TRUE, FALSE))
})

test_that("^ and $ can be used simultaneously", {
  expect_equal(sq_dna %has% "^AGCTAGAACTCTTGATGG$",
               c(FALSE, FALSE, FALSE))
  expect_equal(sq_rna %has% "^GGCUGCGGGACUGAGGC$",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_atp %has% "^mAmYmY$",
               c(TRUE, FALSE, FALSE))
})

# INTERPRETATION AS SUM OF MOTIFS ----
test_that("%has% of many motifs is equal to logical AND of many %has% with one motif", {
  expect_equal(
    sq_ami %has% c("GI", "MATE"),
    (sq_ami %has% "GI") & (sq_ami %has% "MATE")
  )
  expect_equal(
    sq_dna %has% c("CC", "TAAT", "ATGC"),
    (sq_dna %has% "CC") & (sq_dna %has% "TAAT") & (sq_dna %has% "ATGC")
  )
  expect_equal(
    sq_rna %has% c("GC$", "GA"),
    (sq_rna %has% "GC$") & (sq_rna %has% "GA")
  )
  expect_equal(
    sq_unt %has% c("@", "%"),
    (sq_unt %has% "@") & (sq_unt %has% "%")
  )
  expect_equal(
    sq_atp %has% c("mY", "nsA"),
    (sq_atp %has% "mY") & (sq_atp %has% "nsA")
  )
})
