check_sq_type <- function(type, null.ok = FALSE, unt.ok = FALSE, atp.ok = FALSE) {
  # TODO: rethink the idea
  check_choice(type,
               choices = c("dna_bsc", "dna_ext", "rna_bsc",
                           "rna_ext", "ami_bsc", "ami_ext",
                           if (unt.ok) "unt",
                           if (atp.ok) "atp"),
               null.ok = null.ok)
}
assert_sq_type <- makeAssertionFunction(check_sq_type)
test_sq_type <- makeTestFunction(check_sq_type)
expect_sq_type <- makeExpectationFunction(check_sq_type)

check_package_installed <- function(package) {
  if (!package %in% rownames(installed.packages()))
    paste0("you need to install '", package, "' package")
  else
    TRUE
}
assert_package_installed <- makeAssertionFunction(check_package_installed)
test_package_installed <- makeTestFunction(check_package_installed)
expect_package_installed <- makeExpectationFunction(check_package_installed)

check_warning_handling <- function(method) {
  check_choice(method,
               choices = c("error", "warning", "message", "none"))
}
assert_warning_handling <- makeAssertionFunction(check_warning_handling)
test_warning_handling <- makeTestFunction(check_warning_handling)
expect_warning_handling <- makeExpectationFunction(check_warning_handling)

check_motifs_for_type <- function(motifs, type) {
  check_subset(unlist(strsplit(motifs, "")),
               choices = c(get_standard_alphabet(type), "^", "$"))
}
assert_motifs_for_type <- makeAssertionFunction(check_motifs_for_type)
test_motifs_for_type <- makeTestFunction(check_motifs_for_type)
expect_motifs_for_type <- makeExpectationFunction(check_motifs_for_type)

check_alph_regex_friendly <- function(alph) {
  # TODO: reconsider name, it's pretty poor regex now
  check_disjunct(alph, c("^", "$"))
}
assert_alph_regex_friendly <- makeAssertionFunction(check_alph_regex_friendly)
test_alph_regex_friendly <- makeTestFunction(check_alph_regex_friendly)
expect_alph_regex_friendly <- makeExpectationFunction(check_alph_regex_friendly)
