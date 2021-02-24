# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_rna <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")

sq_dna_bsc <- sq(str_dna, alphabet = "dna_bsc")
sq_dna_ext <- sq(str_dna, alphabet = "dna_ext")
sq_rna <- sq(str_rna, alphabet = "rna_ext")
sq_ami <- sq(str_ami, alphabet = "ami_ext")

names_dna <- c("vengeance", "is", "never", "a", "rubber_duck")
names_rna <- c("Monza", "Imola", "Mugello", "Pescara", "Modena")
names_ami <- c("proteins", "vitamins", "fats")

bioseq_dna <- bioseq::new_dna(str_dna)
bioseq_dna_n <- bioseq::new_dna(setNames(str_dna, names_dna))
bioseq_rna <- bioseq::new_rna(str_rna)
bioseq_rna_n <- bioseq::new_rna(setNames(str_rna, names_rna))
bioseq_ami <- bioseq::new_aa(str_ami)
bioseq_ami_n <- bioseq::new_aa(setNames(str_ami, names_ami))

# IMPORT ----
test_that("correctly imports bioseq::bioseq_dna", {
  expect_identical(import_sq(bioseq_dna)[["sq"]],
                   sq_dna_ext)
  expect_identical(import_sq(bioseq_dna_n)[["sq"]],
                   sq_dna_ext)
  expect_identical(import_sq(bioseq_dna_n)[["name"]],
                   names_dna)
})

test_that("correctly imports bioseq::bioseq_rna", {
  expect_identical(import_sq(bioseq_rna)[["sq"]],
                   sq_rna)
  expect_identical(import_sq(bioseq_rna_n)[["sq"]],
                   sq_rna)
  expect_identical(import_sq(bioseq_rna_n)[["name"]],
                   names_rna)
})

test_that("correctly imports bioseq::bioseq_aa", {
  expect_identical(import_sq(bioseq_ami)[["sq"]],
                   sq_ami)
  expect_identical(import_sq(bioseq_ami_n)[["sq"]],
                   sq_ami)
  expect_identical(import_sq(bioseq_ami_n)[["name"]],
                   names_ami)
})

# EXPORT ----
test_that("correctly exports sq object to bioseq::bioseq_dna", {
  expect_identical(export_sq(sq_dna_bsc, "bioseq::bioseq_dna"),
                   bioseq_dna)
  expect_identical(export_sq(sq_dna_ext, "bioseq::bioseq_dna"),
                   bioseq_dna)
  expect_identical(export_sq(sq_dna_bsc, "bioseq::bioseq_dna", name = names_dna),
                   bioseq_dna_n)
})

test_that("correctly exports sq object to bioseq::bioseq_rna", {
  expect_identical(export_sq(sq_rna, "bioseq::bioseq_rna"),
                   bioseq_rna)
  expect_identical(export_sq(sq_rna, "bioseq::bioseq_rna", name = names_rna),
                   bioseq_rna_n)
})

test_that("correctly exports sq object to bioseq::bioseq_aa", {
  expect_identical(export_sq(sq_ami, "bioseq::bioseq_aa"),
                   bioseq_ami)
  expect_identical(export_sq(sq_ami, "bioseq::bioseq_aa", name = names_ami),
                   bioseq_ami_n)
})
