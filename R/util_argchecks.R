.check_simple <- function(check, argname, msg) {
  if (check) stop(argname, " ", msg, call. = FALSE)
}

.check_all_elem <- function(check, argname, msg) {
  if (check) stop("all elements of ", argname, " ", msg, call. = FALSE)
}



.check_isnt_missing <- function(obj, argname) {
  .check_simple(missing(obj), argname, "is missing")
}

.check_isnt_null <- function(obj, argname) {
  .check_simple(is.null(obj), argname, "cannot be NULL")
}

.check_is_single_elem <- function(obj, argname) {
  .check_simple(length(obj) != 1, argname, "should have length 1")
}

.check_isnt_zero_len <- function(obj, argname) {
  .check_simple(length(obj) == 0, argname, "cannot have length 0")
}

.check_has_no_na <- function(obj, argname) {
  .check_simple(any(is.na(obj)), argname, "cannot contain NA")
}

.check_is_unique <- function(obj, argname) {
  .check_simple(length(unique(obj)) != length(obj), argname, "has to have unique elements")
}

.check_is_named <- function(obj, argname) {
  .check_simple(is.null(names(obj)), argname, "should be named")
}

.check_class_character <- function(obj, argname) {
  .check_simple(!is.character(obj), argname, "has to be character")
}

.check_class_logical <- function(obj, argname) {
  .check_simple(!is.logical(obj), argname, "has to be logical")
}

.check_class_numeric <- function(obj, argname) {
  .check_simple(!is.numeric(obj), argname, "has to be numeric")
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

.check_character <- function(obj, argname, ...) {
  .standard_checks(obj, argname, ...)
  if (!is.null(obj)) .check_class_character(obj, argname)
}

.check_logical <- function(obj, argname, ...) {
  .standard_checks(obj, argname, ...)
  if (!is.null(obj)) .check_class_logical(obj, argname)
}

.check_numeric <- function(obj, argname, ..., 
                           allow_nan = FALSE, allow_inf = FALSE, allow_negative = FALSE,
                           allow_zero = TRUE) {
  .standard_checks(obj, argname, ...)
  if (!is.null(obj)) {
    .check_class_numeric(obj, argname)
    if (!allow_nan )     .check_simple(any(is.nan(obj)),       argname, "cannot contain NaN values")
    if (!allow_inf )     .check_simple(any(is.infinite(obj)),  argname, "cannot contain infinite values")
    if (!allow_zero)     .check_simple(any(obj == 0),          argname, "cannot be equal to 0")
    if (!allow_negative) .check_simple(any(obj < 0),           argname, "cannot be negative")
  }
}

.check_integer <- function(obj, argname, ...) {
  .check_numeric(obj, argname, ...)
  if (!is.null(obj)) 
    .check_simple(any(floor(obj) != obj), argname, "has to be integer")
}

.check_type <- function(obj, argname = "'type'", allow_null = FALSE, allow_unt = FALSE) {
  .check_isnt_missing(obj, argname)
  if (!allow_null) .check_isnt_null(obj, argname)
  else if (!is.null(obj)) {
    allowed <- c("ami", "nuc", if (allow_unt) "unt")
    .check_simple(!obj %in% allowed, argname, paste0("has to be one of '", paste(allowed, collapse = "', '"), "'")) 
  }
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

.check_eq_lens <- function(obj_1, obj_2, name_1, name_2){
  .check_simple(length(obj_1) != length(obj_2), paste0(name_1, " and ", name_2), "have to have equal lengths")
}



### specific checks - used mainly once

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

.check_export_format <- function(export_format, ami_formats, nuc_formats) {
  if (missing(export_format) ||
      !(export_format %in% c(ami_formats, nuc_formats))) {
    stop("you need to specify proper 'export_format'; check manual for possible formats", call. = FALSE)
  }
}

.check_type_matches_format <- function(type, export_format, ami_formats, nuc_formats) {
  if (type == "ami" && !(export_format %in% ami_formats) ||
      type == "nuc" && !(export_format %in% nuc_formats)) 
    stop("'sq' object type doesn't match 'export_format'", call. = FALSE)
}

.check_enc_names_in_alph <- function(encoding, alph) {
  if (!all(names(encoding) %in% alph)) 
    stop("all names of 'encoding' has to be letters from alphabet (elements of 'alphabet' attribute of 'sq')", call. = FALSE)
}
  
.check_has_both_UT <- function(has_U, has_T) {
  if (has_U && has_T)
    stop("'nucsq' sequences contains both 'U' and 'T' letters - should contain only one of them", call. = FALSE)
}

.check_has_A_no_UT <- function(has_U, has_T, has_A) {
  if (!has_U && !has_T && has_A)
    stop("'nucsq' sequences contains 'A' elements, but does not contain 'T' nor 'U' - unable to guess if it's dna or rna", call. = FALSE)
}

.check_is_dna_matches <- function(is_dna, has_U, has_T) {
  if ((is_dna && has_U) || (!is_dna && has_T)) 
    stop("if 'is_dna' is TRUE, sequences cannot contain 'U'; if is FALSE, sequences cannot contain 'T'", call. = FALSE)
}

.check_motifs_proper_alph <- function(motifs, type, alph = NULL) {
  if (type %in% c("ami", "nuc")) {
    if (!all(unlist(strsplit(motifs, "")) %in% c(.get_standard_alph(type, FALSE), "^", "$"))) 
      stop("motifs that you're searching for in the 'sq' object needs to consist of letters from its alphabet and optionally '^' or '$' characters", call. = FALSE)
  } else if (any(alph %in% c("^", "$", "?", "(", "=", ")", "\\", ".", "|", "+", "*", "{", "}", "[", "]"))) 
    stop("you cannot search for motifs if any of those characters: ^$?=()\\.|+*{}[] are elements of 'sq' alphabet; if you want to use them, please substitute those letters with some other using 'substitute_letters'", call. = FALSE)
}

.check_all_up_alph_proper <- function(up_alph, dest_alph) {
  if (!all(up_alph %in% dest_alph)) 
    stop("some sequences have levels that are invalid for given 'dest_type'; you can check them with 'get_invalid_letters' function and fix them with 'substitute_letters'", call. = FALSE)
}

.check_paste_or_na <- function(paste_char, use_na_char) {
  if (paste_char && !use_na_char) 
    stop("'paste_char' can be TRUE if and only if 'use_na_char' is FALSE", call. = FALSE)
}

.check_list_dists <- function(dists) {
  .check_simple(!all(sapply(dists, is.numeric)), "each element of 'dists' has to be integer")
  dists_f <- unlist(dists)
  .check_isnt_null(dists_f, "any of elements of 'dists'")
  .check_has_no_na(dists_f, "any of elements of 'dists'")
  .check_simple(any(is.nan(dists_f)), "any of elements of 'dists'", "cannot contain NaN values")
  .check_simple(any(is.infinite(dists_f)), "any of elements of 'dists'", "cannot contain infinite values")
  .check_simple(any(dists_f < 0), "any of elements of 'dists'", "cannot be negative")
  #.check_is_unique(dists, "'dists'")
}

.check_alph_is_subset <- function(sq, alph) {
  if (!all(alph %in% .get_alph(sq)))
    stop("'alph' contains letters that aren't elements of alphabet of 'sq'", call. = FALSE)
}

.check_dists_prop_len <- function(sq, dists) {
  min_sq_len <- min(.get_lens(sq))
  if (is.list(dists)) {
    if (any(sapply(dists, function(d) sum(d) + length(d) + 1) > min_sq_len))
      stop("some sequences in 'sq' are shorter than some of kmers to extract ", call. = FALSE)
  } else if (sum(dists) + length(dists) + 1 > min_sq_len)
    stop("some sequences in 'sq' are shorter than given kmers", call. = FALSE)
}