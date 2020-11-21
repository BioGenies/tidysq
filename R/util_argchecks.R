assert_sq_type <- function(type, null.ok = FALSE, unt.ok = FALSE) {
  # TODO: rethink the idea
  assert_choice(type,
                choices = c("dna_bsc", "dna_ext", "rna_bsc", "rna_ext",
                            "ami_bsc", "ami_ext", if (unt.ok) "unt"),
                null.ok = null.ok)
}

assert_package_installed <- function(package) {
  if (!package %in% rownames(installed.packages()))
    stop("you need to install '", package, "' package to export object to its formats", call. = FALSE)
  invisible(package)
}

assert_motifs_for_type <- function(motifs, type) {
  assert_subset(unlist(strsplit(motifs, "")), c(CPP_get_standard_alphabet(type), "^", "$"))
}

assert_alph_regex_friendly <- function(alph) {
  # TODO: reconsider name, it's pretty poor regex now
  assert_disjunct(alph, c("^", "$"))
}
