# general checks - used basically everywhere ----
assert_sq_type <- function(type, null.ok = FALSE, unt.ok = FALSE) {
  assert_choice(type,
                choices = c("ami", "dna", "rna", if (unt.ok) "unt"),
                null.ok = null.ok)
}

.check_simple <- function(check, argname, msg) {
  if (check) stop(argname, " ", msg, call. = FALSE)
}

.check_all_elem <- function(check, argname, msg) {
  if (check) stop("all elements of ", argname, " ", msg, call. = FALSE)
}



.standard_checks <- function(obj, argname,
                             allow_null = FALSE, 
                             single_elem = FALSE, 
                             allow_zero_len = FALSE,
                             allow_na = FALSE) {
  
                       .check_isnt_missing(obj, argname)
  if (!allow_null    ) .check_isnt_null(obj, argname) else if (is.null(obj)) return()
  if (single_elem    ) .check_is_single_elem(obj, argname)
  if (!allow_zero_len) .check_isnt_zero_len(obj, argname)
  if (!allow_na      ) .check_has_no_na(obj, argname)
}

.check_sq_has_type <- function(sq, argname, type) {
  .check_simple(!type %in% .get_sq_type(sq), argname, paste0("has to have '", type, "' type"))
}

.check_is_clean <- function(sq, argname) {
  .check_simple(!.is_cleaned(sq), argname, "has to be clean (has to have 'cln' subtype)")
}


.check_nchar <- function(obj, argname, allow_zero_nchar = FALSE, requested_nchar = NULL,
                         demand_eq_len = FALSE, minimal_nchar = NULL) {
  if (!is.null(requested_nchar)) 
    .check_all_elem(!all(nchar(obj) == requested_nchar), argname, 
                    paste0("have to have length of ", requested_nchar))
  else if (!is.null(minimal_nchar)) 
    .check_all_elem(!all(nchar(obj) >= requested_nchar), argname, 
                    paste0("have to have length at least equal to ", minimal_nchar))
  if (!allow_zero_nchar) 
    .check_all_elem(!all(nchar(obj) != 0), argname, "have to have positive length")
  if (demand_eq_len) {
    lens <- nchar()
    .check_all_elem(length(unique(lens)) == 1, argname, "have to have equal length")
  }
}

.check_type_or_nonst_alph <- function(type, is_clean, non_standard) {
  if (!(is.null(type) && 
        is.null(is_clean)) &&
      !is.null(non_standard))
    stop("if you specify 'non_standard', you cannot specify neither 'type' nor 'is_clean'", call. = FALSE)
}

.check_alph_length <- function(alph) {
  if (length(alph) > 30) 
    stop("max length of alphabet is 30 letters, sequences that are being constructed exceed this limit", call. = FALSE)
}


# specific checks - used mainly once ----

.check_real_alph_clean <- function(real_alph, type, is_clean) {
  if (!is.null(is_clean) &&
      is_clean == TRUE &&
      !all(real_alph %in% .get_standard_alph(type, TRUE))) {
    stop("'is_clean' is given TRUE, but sequences contain at least one ambiguous element", call. = FALSE)
  }
}

.check_alph_matches_type <- function(alph, type, is_clean = NULL) {
  if (!is.null(type) && !(type == "unt")) {
    if (is.null(is_clean)) is_clean <- FALSE
    if (!all(alph %in% .get_standard_alph(type, is_clean))) {
      stop("there are letters in alphabet that aren't suitable for given 'type' or 'is_clean' parameters", call. = FALSE)
    }
  }
}

.check_is_installed <- function(package) {
  if (!package %in% rownames(installed.packages()))
      stop("you need to install '", package, "' package to export object to its formats", call. = FALSE)
}

.check_motifs_proper_alph <- function(motifs, type, alph = NULL) {
  if (type %in% c("ami", "dna", "rna")) {
    if (!all(unlist(strsplit(motifs, "")) %in% c(.get_standard_alph(type, FALSE), "^", "$"))) 
      stop("motifs that you're searching for in the 'sq' object needs to consist of letters from its alphabet and optionally '^' or '$' characters", call. = FALSE)
  } else if (any(alph %in% c("^", "$", "?", "(", "=", ")", "\\", ".", "|", "+", "*", "{", "}", "[", "]"))) 
    stop("you cannot search for motifs if any of those characters: ^$?=()\\.|+*{}[] are elements of 'sq' alphabet; if you use them, please substitute those letters with some other using 'substitute_letters'", call. = FALSE)
}
