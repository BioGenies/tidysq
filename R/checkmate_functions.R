check_sq_type <- function(type, null.ok = FALSE, unt.ok = FALSE) {
  # TODO: rethink the idea
  check_choice(type,
               choices = c("dna_bsc", "dna_ext", "rna_bsc", "rna_ext",
                           "ami_bsc", "ami_ext", if (unt.ok) "unt"),
               null.ok = null.ok)
}

assert_sq_type <- makeAssertionFunction(check_sq_type)

expect_sq_type <- makeExpectationFunction(check_sq_type)

check_package_installed <- function(package) {
  if (!package %in% rownames(installed.packages()))
    paste0("you need to install '", package, "' package")
  else
    TRUE
}

assert_package_installed <- makeAssertionFunction(check_package_installed)

expect_package_installed <- makeExpectationFunction(check_package_installed)

assert_motifs_for_type <- function(motifs, type) {
  assert_subset(unlist(strsplit(motifs, "")), c(get_standard_alphabet(type), "^", "$"))
}

assert_alph_regex_friendly <- function(alph) {
  # TODO: reconsider name, it's pretty poor regex now
  assert_disjunct(alph, c("^", "$"))
}
