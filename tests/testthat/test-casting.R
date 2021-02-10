# SETUP ----
str_dna_bsc <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_dna_ext <- c("NARTYVTCY", "", "ATKCYGDD", "", "DNAKYTD")
str_rna_bsc <- c("AUCAGU-A-GU-CA", "CUG-A-CUGAG-CC", "-CUG-AGAGU-UA-")
str_rna_ext <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_ami_bsc <- c("ACEH", "PASAI", "MALACCA", "SIAK")
str_ami_ext <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
str_unt <- c("VIP01", "VIP02", "VIP04", "MISSING_ONE")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

matrix_rna_bsc <- matrix(
  c("A", "U", "C", "A", "G", "U", "-", "A", "-", "G", "U", "-", "C", "A",
    "C", "U", "G", "-", "A", "-", "C", "U", "G", "A", "G", "-", "C", "C",
    "-", "C", "U", "G", "-", "A", "G", "A", "G", "U", "-", "U", "A", "-"),
  nrow = 3, byrow = TRUE
)
matrix_ami_bsc <- matrix(
  c("A", "C", "E", "H", NA, NA, NA,
    "P", "A", "S", "A", "I", NA, NA,
    "M", "A", "L", "A", "C", "C", "A",
    "S", "I", "A", "K", NA, NA, NA),
  nrow = 4, byrow = TRUE
)
matrix_atp <- matrix(
  c("mA", "mY", "mY",
    "nbA", "nsA", "mA",
    NA, NA, NA),
  nrow = 3, byrow = TRUE
)

alph_atp <- c("mA", "mY", "nbA", "nsA")

sq_dna_bsc <- sq(str_dna_bsc, alphabet = "dna_bsc")
sq_dna_ext <- sq(str_dna_ext, alphabet = "dna_ext")
sq_rna_bsc <- sq(str_rna_bsc, alphabet = "rna_bsc")
sq_rna_ext <- sq(str_rna_ext, alphabet = "rna_ext")
sq_ami_bsc <- sq(str_ami_bsc, alphabet = "ami_bsc")
sq_ami_ext <- sq(str_ami_ext, alphabet = "ami_ext")
sq_unt <- sq(str_unt, alphabet = "unt")
sq_atp <- sq(str_atp, alphabet = alph_atp)

sq_with_na <- sq(str_dna_ext, alphabet = "dna_bsc")
str_with_na_1 <- c("!A!T!!TC!", "", "AT!C!G!!", "", "!!A!!T!")
str_with_na_2 <- c("?A?T??TC?", "", "AT?C?G??", "", "??A??T?")
str_with_na_3 <- c("<?>A<?>T<?><?>TC<?>", "", "AT<?>C<?>G<?><?>", "", "<?><?>A<?><?>T<?>")

# added because of Biostrings warning
suppressWarnings({
  biostr_dna_bsc <- Biostrings::DNAStringSet(str_dna_bsc)
})

seqinr_ami_bsc <- lapply(str_ami_bsc, function(x)
  seqinr::as.SeqFastaAA(seqinr::s2c(x)))

# CASTING TO SQ ----
test_that("character vector is casted to sq with as.sq()", {
  expect_identical(as.sq(str_rna_bsc),
                   sq(str_rna_bsc))
  expect_identical(as.sq(str_unt),
                   sq(str_unt))
  expect_identical(as.sq(str_atp),
                   sq(str_atp))
})

test_that("non-character objects are passed to import_sq()", {
  expect_identical(as.sq(biostr_dna_bsc),
                   import_sq(biostr_dna_bsc))
  expect_identical(as.sq(seqinr_ami_bsc),
                   import_sq(seqinr_ami_bsc))
  expect_error(as.sq(function(x, y) x + y))
})

# CASTING TO CHARACTER ----
test_that("applying as.character() returns original character vector", {
  expect_equivalent(as.character(sq_dna_bsc), str_dna_bsc)
  expect_equivalent(as.character(sq_dna_ext), str_dna_ext)
  expect_equivalent(as.character(sq_rna_bsc), str_rna_bsc)
  expect_equivalent(as.character(sq_rna_ext), str_rna_ext)
  expect_equivalent(as.character(sq_ami_bsc), str_ami_bsc)
  expect_equivalent(as.character(sq_ami_ext), str_ami_ext)
  expect_equivalent(as.character(sq_unt), str_unt)
  expect_equivalent(as.character(sq_atp), str_atp)
})

test_that("as.character() has NA_letter parameter with default value", {
  withr::local_options(list(tidysq_NA_letter = "!"))
  expect_equivalent(as.character(sq_with_na), str_with_na_1)
  expect_equivalent(as.character(sq_with_na, NA_letter = "?"), str_with_na_2)
  expect_equivalent(as.character(sq_with_na, NA_letter = "<?>"), str_with_na_3)
})

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

# CASTING TO MATRIX ----
test_that("as.matrix() creates a character matrix with a row for each sequence and a column for each element", {
  expect_matrix(as.matrix(sq_rna_bsc),
                mode = "character",
                nrows = vec_size(sq_rna_bsc),
                ncols = max(get_sq_lengths(sq_rna_bsc)))
  expect_matrix(as.matrix(sq_dna_ext),
                mode = "character",
                nrows = vec_size(sq_dna_ext),
                ncols = max(get_sq_lengths(sq_dna_ext)))
  expect_matrix(as.matrix(sq_unt),
                mode = "character",
                nrows = vec_size(sq_unt),
                ncols = max(get_sq_lengths(sq_unt)))
  expect_matrix(as.matrix(sq_atp),
                mode = "character",
                nrows = vec_size(sq_atp),
                ncols = max(get_sq_lengths(sq_atp)))
})

test_that("as.matrix() splits equal-length sequences into columns", {
  expect_identical(as.matrix(sq_rna_bsc),
                   matrix_rna_bsc)
})

test_that("as.matrix() fills shorter sequences with NAs at the end", {
  expect_identical(as.matrix(sq_ami_bsc),
                   matrix_ami_bsc)
})

test_that("as.matrix() splits sequences by letter, not by character", {
  expect_identical(as.matrix(sq_atp),
                   matrix_atp)
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
