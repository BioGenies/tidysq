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

.check_motifs_proper_alph <- function(motifs, type, alph = NULL) {
  if (type %in% c("ami", "dna", "rna")) {
    if (!all(unlist(strsplit(motifs, "")) %in% c(get_standard_alphabet(type), "^", "$")))
      stop("motifs that you're searching for in the 'sq' object needs to consist of letters from its alphabet and optionally '^' or '$' characters", call. = FALSE)
  } else if (any(alph %in% c("^", "$", "?", "(", "=", ")", "\\", ".", "|", "+", "*", "{", "}", "[", "]"))) 
    stop("you cannot search for motifs if any of those characters: ^$?=()\\.|+*{}[] are elements of 'sq' alphabet; if you use them, please substitute those letters with some other using 'substitute_letters'", call. = FALSE)
}
