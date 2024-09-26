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

# CONCATENATION OF SAME TYPE SEQUENCES ----
test_that("c() on same type sequences return object of the same prototype", {
  expect_vector(c(sq_dna_bsc, sq_dna_bsc),
                ptype = vec_ptype(sq_dna_bsc),
                size = 2 * vec_size(sq_dna_bsc))
  expect_vector(c(sq_dna_ext, sq_dna_ext),
                ptype = vec_ptype(sq_dna_ext),
                size = 2 * vec_size(sq_dna_ext))
  expect_vector(c(sq_rna_bsc, sq_rna_bsc),
                ptype = vec_ptype(sq_rna_bsc),
                size = 2 * vec_size(sq_rna_bsc))
  expect_vector(c(sq_rna_ext, sq_rna_ext),
                ptype = vec_ptype(sq_rna_ext),
                size = 2 * vec_size(sq_rna_ext))
  expect_vector(c(sq_ami_bsc, sq_ami_bsc),
                ptype = vec_ptype(sq_ami_bsc),
                size = 2 * vec_size(sq_ami_bsc))
  expect_vector(c(sq_ami_ext, sq_ami_ext),
                ptype = vec_ptype(sq_ami_ext),
                size = 2 * vec_size(sq_ami_ext))
  expect_vector(c(sq_unt, sq_unt),
                ptype = vec_ptype(sq_unt),
                size = 2 * vec_size(sq_unt))
  expect_vector(c(sq_atp, sq_atp),
                ptype = vec_ptype(sq_atp),
                size = 2 * vec_size(sq_atp))
})

test_that("c() on same type sequences return correct value", {
  expect_identical(c(sq_dna_bsc, sq_dna_bsc),
                   sq(c(str_dna_bsc, str_dna_bsc), alphabet = "dna_bsc"))
  expect_identical(c(sq_dna_ext, sq_dna_ext),
                   sq(c(str_dna_ext, str_dna_ext), alphabet = "dna_ext"))
  expect_identical(c(sq_rna_bsc, sq_rna_bsc),
                   sq(c(str_rna_bsc, str_rna_bsc), alphabet = "rna_bsc"))
  expect_identical(c(sq_rna_ext, sq_rna_ext),
                   sq(c(str_rna_ext, str_rna_ext), alphabet = "rna_ext"))
  expect_identical(c(sq_ami_bsc, sq_ami_bsc),
                   sq(c(str_ami_bsc, str_ami_bsc), alphabet = "ami_bsc"))
  expect_identical(c(sq_ami_ext, sq_ami_ext),
                   sq(c(str_ami_ext, str_ami_ext), alphabet = "ami_ext"))
  expect_identical(c(sq_unt, sq_unt),
                   sq(c(str_unt, str_unt), alphabet = "unt"))
  expect_identical(c(sq_atp, sq_atp),
                   sq(c(str_atp, str_atp), alphabet = alph_atp))
})

# CONCATENATION OF BASIC AND EXTENDED TYPES ----
test_that("c() generalizes basic and extended types to extended", {
  expect_vector(c(sq_dna_bsc, sq_dna_ext, sq_dna_bsc),
                ptype = vec_ptype(sq_dna_ext),
                size = 2 * vec_size(sq_dna_bsc) + vec_size(sq_dna_ext))
  expect_vector(c(sq_rna_bsc, sq_rna_ext, sq_rna_bsc),
                ptype = vec_ptype(sq_rna_ext),
                size = 2 * vec_size(sq_rna_bsc) + vec_size(sq_rna_ext))
  expect_vector(c(sq_ami_bsc, sq_ami_ext, sq_ami_bsc),
                ptype = vec_ptype(sq_ami_ext),
                size = 2 * vec_size(sq_ami_bsc) + vec_size(sq_ami_ext))
})

test_that("c() on basic and extended sequences return correct value", {
  expect_identical(c(sq_dna_bsc, sq_dna_ext, sq_dna_bsc),
                   sq(c(str_dna_bsc, str_dna_ext, str_dna_bsc), alphabet = "dna_ext"))
  expect_identical(c(sq_rna_bsc, sq_rna_ext, sq_rna_bsc),
                   sq(c(str_rna_bsc, str_rna_ext, str_rna_bsc), alphabet = "rna_ext"))
  expect_identical(c(sq_ami_bsc, sq_ami_ext, sq_ami_bsc),
                   sq(c(str_ami_bsc, str_ami_ext, str_ami_bsc), alphabet = "ami_ext"))
})

# CONCATENATION WITH UNT SQ ----
test_that("c() generalizes any standard type and untyped sequences to unt sq", {
  expect_vector(c(sq_dna_bsc, sq_unt, sq_dna_bsc),
                ptype = sq_ptype(union(as.character(alphabet(sq_dna_bsc)), as.character(alphabet(sq_unt))), "unt"),
                size = 2 * vec_size(sq_dna_bsc) + vec_size(sq_unt))
  expect_vector(c(sq_dna_ext, sq_unt, sq_dna_ext),
                ptype = sq_ptype(union(as.character(alphabet(sq_dna_ext)), as.character(alphabet(sq_unt))), "unt"),
                size = 2 * vec_size(sq_dna_ext) + vec_size(sq_unt))
  expect_vector(c(sq_rna_bsc, sq_unt, sq_rna_bsc),
                ptype = sq_ptype(union(as.character(alphabet(sq_rna_bsc)), as.character(alphabet(sq_unt))), "unt"),
                size = 2 * vec_size(sq_rna_bsc) + vec_size(sq_unt))
  expect_vector(c(sq_rna_ext, sq_unt, sq_rna_ext),
                ptype = sq_ptype(union(as.character(alphabet(sq_rna_ext)), as.character(alphabet(sq_unt))), "unt"),
                size = 2 * vec_size(sq_rna_ext) + vec_size(sq_unt))
  expect_vector(c(sq_ami_bsc, sq_unt, sq_ami_bsc),
                ptype = sq_ptype(union(as.character(alphabet(sq_ami_bsc)), as.character(alphabet(sq_unt))), "unt"),
                size = 2 * vec_size(sq_ami_bsc) + vec_size(sq_unt))
  expect_vector(c(sq_ami_ext, sq_unt, sq_ami_ext),
                ptype = sq_ptype(union(as.character(alphabet(sq_ami_ext)), as.character(alphabet(sq_unt))), "unt"),
                size = 2 * vec_size(sq_ami_ext) + vec_size(sq_unt))
})

test_that("c() on any standard type and untyped sequences return correct value", {
  expect_identical(c(sq_dna_bsc, sq_unt, sq_dna_bsc),
                   {
                     x <- sq(c(str_dna_bsc, str_unt, str_dna_bsc), alphabet = union(as.character(alphabet(sq_dna_bsc)), as.character(alphabet(sq_unt))))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
  expect_identical(c(sq_dna_ext, sq_unt, sq_dna_ext),
                   {
                     x <- sq(c(str_dna_ext, str_unt, str_dna_ext), alphabet = union(as.character(alphabet(sq_dna_ext)), as.character(alphabet(sq_unt))))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
  expect_identical(c(sq_rna_bsc, sq_unt, sq_rna_bsc),
                   {
                     x <- sq(c(str_rna_bsc, str_unt, str_rna_bsc), alphabet = union(as.character(alphabet(sq_rna_bsc)), as.character(alphabet(sq_unt))))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
  expect_identical(c(sq_rna_ext, sq_unt, sq_rna_ext),
                   {
                     x <- sq(c(str_rna_ext, str_unt, str_rna_ext), alphabet = union(as.character(alphabet(sq_rna_ext)), as.character(alphabet(sq_unt))))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
  expect_identical(c(sq_ami_bsc, sq_unt, sq_ami_bsc),
                   {
                     x <- sq(c(str_ami_bsc, str_unt, str_ami_bsc), alphabet = union(as.character(alphabet(sq_ami_bsc)), as.character(alphabet(sq_unt))))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
  expect_identical(c(sq_ami_ext, sq_unt, sq_ami_ext),
                   {
                     x <- sq(c(str_ami_ext, str_unt, str_ami_ext), alphabet = union(as.character(alphabet(sq_ami_ext)), as.character(alphabet(sq_unt))))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
})

# CONCATENATION WITH CHARACTER VECTOR ----
test_that("c() coerces strings to sq if after an sq object and keeps original sq type", {
  expect_vector(c(sq_dna_ext, str_dna_bsc),
                ptype = vec_ptype(sq_dna_ext),
                size = vec_size(sq_dna_ext) + vec_size(str_dna_bsc))
  expect_vector(c(sq_dna_bsc, str_dna_ext),
                ptype = vec_ptype(sq_dna_bsc),
                size = vec_size(sq_dna_bsc) + vec_size(str_dna_ext))
  expect_vector(c(sq_rna_ext, str_rna_bsc),
                ptype = vec_ptype(sq_rna_ext),
                size = vec_size(sq_rna_ext) + vec_size(str_rna_bsc))
  expect_vector(c(sq_rna_bsc, str_rna_ext),
                ptype = vec_ptype(sq_rna_bsc),
                size = vec_size(sq_rna_bsc) + vec_size(str_rna_ext))
  expect_vector(c(sq_ami_ext, str_ami_bsc),
                ptype = vec_ptype(sq_ami_ext),
                size = vec_size(sq_ami_ext) + vec_size(str_ami_bsc))
  expect_vector(c(sq_ami_bsc, str_ami_ext),
                ptype = vec_ptype(sq_ami_bsc),
                size = vec_size(sq_ami_bsc) + vec_size(str_ami_ext))
  expect_vector(c(sq_unt, str_rna_ext),
                ptype = vec_ptype(sq_unt),
                size = vec_size(sq_unt) + vec_size(str_rna_ext))
  expect_vector(c(sq_atp, str_atp),
                ptype = vec_ptype(sq_atp),
                size = vec_size(sq_atp) + vec_size(str_atp))
})

test_that("c() with sq and character vector return correct value", {
  expect_identical(c(sq_dna_ext, str_dna_bsc),
                   sq(c(str_dna_ext, str_dna_bsc), alphabet = "dna_ext"))
  expect_identical(c(sq_dna_bsc, str_dna_ext),
                   sq(c(str_dna_bsc, str_dna_ext), alphabet = "dna_bsc"))
  expect_identical(c(sq_rna_ext, str_rna_bsc),
                   sq(c(str_rna_ext, str_rna_bsc), alphabet = "rna_ext"))
  expect_identical(c(sq_rna_bsc, str_rna_ext),
                   sq(c(str_rna_bsc, str_rna_ext), alphabet = "rna_bsc"))
  expect_identical(c(sq_ami_ext, str_ami_bsc),
                   sq(c(str_ami_ext, str_ami_bsc), alphabet = "ami_ext"))
  expect_identical(c(sq_ami_bsc, str_ami_ext),
                   sq(c(str_ami_bsc, str_ami_ext), alphabet = "ami_bsc"))
  expect_identical(c(sq_unt, str_rna_ext),
                   {
                     x <- sq(c(str_unt, str_rna_ext), alphabet = as.character(alphabet(sq_unt)))
                     attr(alphabet(x), "type") <- "unt"
                     class(x)[class(x) == "sq_atp"] <- "sq_unt"
                     x
                   })
  expect_identical(c(sq_atp, str_atp),
                   sq(c(str_atp, str_atp), alphabet = alph_atp))
})
