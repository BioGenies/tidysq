.check_type_in_ami_nuc_unt_NULL <- function(type) {
  if (!is.null(type) &&
      !type %in% c("unt", "ami", "nuc")) 
    stop("'type' has to be NULL or one of 'unt', 'ami', 'nuc'")
}

.check_type_in_ami_nuc_unt <- function(type) {
  if (!type %in% c("unt", "ami", "nuc")) 
    stop("'type' has to be one of 'unt', 'ami', 'nuc'")
}

.check_nc_type_in_ami_nuc <- function(type) {
  if (missing(type) ||
      !type %in% c("nuc", "ami"))
    stop("in no_check mode 'type' needs to be one of 'nuc', 'ami'")
}

.check_type_in_ami_nuc <- function(type) {
  if (is.null(type) ||
      !type %in% c("nuc", "ami"))
    stop("'type' has to be one of 'nuc', 'ami'")
}

.check_sq_type_in_ami_nuc <- function(type) {
  if (!type %in% c("nuc", "ami"))
    stop("'sq' type has to be one of 'nuc', 'ami'")
}

.check_sqstr_proper_char <- function(sq) {
  if (!is.character(sq) ||
      length(sq) == 0 ||
      any(is.na(sq)) ||
      any(is.null(sq))) 
    stop("'sq' has to be a positive-length vector of non-NA, non-NULL strings")
}

.check_is_clean_in_TRUE_FALSE_NULL <- function(is_clean) {
  if (!is.null(is_clean) &&
      !is_clean %in% c(TRUE, FALSE)) 
    stop("'is_clean' has to be TRUE, FALSE or NULL")
}

.check_nc_is_clean_in_TRUE_FALSE <- function(is_clean) {
  if (!is_clean %in% c(TRUE, FALSE)) 
    stop("'is_clean' has to be TRUE, FALSE in no_ceck mode")
}



.check_nonst_proper_char <- function(non_standard) {
  if (!is.character(non_standard) ||
      length(non_standard) == 0 ||
      any(is.na(non_standard)) ||
      any(is.null(non_standard))) 
    stop("'non_standard' has to be a positive-length character vector without NA")
}

.check_nonst_nchar <- function(non_standard) {
  if (any(nchar(non_standard) < 2)) 
    stop("non standard letters specified in 'non_standard' parameter have to have more than one character")
}



.check_sq_is_clean <- function(sq) {
  if (!.is_cleaned(sq))
    stop("'sq' object has to be clean")
}

.check_inds_are_numeric <- function(indices) {
  if (!(is.numeric(indices) && 
        floor(indices) == indices)) {
    stop("'indices' has to be an integer vector")
  }
}

.check_file_is_char <- function(file) {
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
}

.check_name_proper_char <- function(name) {
  if (missing(name) ||
      is.null(name) ||
      !is.character(name) ||
      any(is.na(name))) {
    stop("'name' has to be a non-NULL character vector without NA's")
  }
}

.check_nchar_proper_int <- function(nchar) {
  if (!is.numeric(nchar) ||
      (floor(nchar) != nchar) ||
      (length(nchar) != 1) ||
      is.na(nchar) || 
      is.nan(nchar) || 
      !is.finite(nchar) ||
      nchar <= 0) 
    stop("'nchar' has to be positive integer indicating max number of elements of sequence in single line in file")
}

.check_eq_lens <- function(sq, name) {
  if (length(sq) != length(name)) 
    stop("'name' has to have length equal to 'sq'")
}

.check_matrix_no_na <- function(matrix) {
  if (any(is.na(matrix)))
    stop("there can't be any NA in 'sq' object")
}

.check_matrix_no_star <- function(matrix) {
  if (any(matrix == "*"))
    stop("'sq' object cannot contain '*'; you can remove them using substitute_letters and remove_na")    
}

.check_num_pos_proper_int <- function(num_pos) {
  if (!is.numeric(num_pos) ||
      (floor(num_pos) != num_pos) ||
      (length(nchar) != 1) ||
      is.na(num_pos) || 
      is.nan(num_pos) || 
      !is.finite(num_pos) ||
      num_pos <= 0) 
    stop("'nchar' has to be positive integer indicating number of rows (number of positions, sequence lenghts) in PSSM")
}

.check_mode_is_proper <- function(mode) {
  if (missing(mode) ||
      !mode %in% c("counts", "freqs", "shannon", "information", "kullback-leibler"))
    stop("'mode' has to be one of 'counts', 'freqs', 'shannon', 'information', 'kullback-leibler'")
}

.check_back_dist_is_proper <- function(background_dist) {
  if (!background_dist %in% rownames(bg_freqs))
    stop("'background_dist' has to be in rownames of BGFREQS dataframe")
}

.check_sq_lens_eq <- function(sq) {
  if (length(unique(.get_lens(sq))) != 1) 
    stop("all sequences in 'sq' have to have the same length")
}

.check_is_pssm <- function(pssm) {
  if (!is.numeric(pssm) || 
      !is.matrix(pssm) ||
      (colnames(pssm) != setdiff(.get_standard_alph("ami", TRUE), "*") &&
       colnames(pssm) != .get_standard_alph("nuc", TRUE)))
    stop("'pssm' has to be pssm matrix (result of compute_pssm function)")
}

.check_is_pssm_or_numeric <- function(pssm) {
  if (!is.numeric(pssm) ||
      (!is.matrix(pssm) ||
       (colnames(pssm) != setdiff(.get_standard_alph("ami", TRUE), "*") &&
        colnames(pssm) != .get_standard_alph("nuc", TRUE))))
    stop("'pssm' has to be pssm matrix (result of compute_pssm function) or numeric vector")
}

.check_types_match <- function(sq, pssm) {
  type <- .get_sq_type(sq)
  alph <- .get_alph(sq)
  if (type == "ami" && !identical(colnames(pssm), setdiff(alph, "*")) ||
      type == "nuc" && !identical(colnames(pssm), alph))
    stop("alphabet of 'sq' object doesn't match 'pssm' alphabet")
}

.check_lens_eq_num_pos <- function(len, num_pos) {
  if (len != num_pos) 
    stop("length of sequences in 'sq' has to be the same as number of rows of 'pssm'")
}

.check_enc_proper_int <- function(encoding) {
  if (length(encoding) == 0 ||
      (floor(encoding) != encoding) ||
      any(is.nan(encoding)) || 
      any(!is.finite(encoding)) ||
      any(encoding < 0)) 
    stop("if 'encoding' is numeric, it has to be non-negative integer vector, with no NaN's or non-finite values")
}

.check_enc_proper_char <- function(encoding) {
  if (is.null(encoding) ||
      length(encoding) == 0 ||
      !is.character(encoding) ||
      any(nchar(encoding) == 0))
    stop("'encoding' has to be either integer or character vector of positive lenght")
}

.check_n_is_int <- function(n) {
  if (!is.numeric(n) ||
      (floor(n) != n) ||
      (length(n) != 1) ||
      is.na(n) || 
      is.nan(n) || 
      !is.finite(n) ||
      n <= 0) 
    stop("'n' has to positive, finite, non-NA, non-NaN integer")
}

.check_len_is_int <- function(len) {
  if (!is.numeric(len) ||
      (floor(len) != len) ||
      (length(len) != 1) ||
      is.na(len) || 
      is.nan(len) || 
      !is.finite(len) ||
      len <= 0) 
    stop("'len' has to positive, finite, non-NA, non-NaN integer")
}
  
.check_is_clean_in_TRUE_FALSE <- function(is_clean) {
  if (is.null(is_clean) ||
      length(is_clean) != 1 ||
      !is_clean %in% c(TRUE, FALSE))
    stop("'is_clean' has to be either TRUE or FALSE")
}

.check_sd_is_numeric_or_NULL <- function(sd) {
  if (!is.null(sd) &&
      (!is.numeric(sd) ||
       (floor(sd) != sd) ||
       (length(sd) != 1) ||
       is.na(sd) || 
       is.nan(sd) || 
       !is.finite(sd) ||
       sd <= 0))
    stop("'sd' has to be NULL or positive, finite, non-NA, non-NaN integer")
}

.check_use_gap_in_TRUE_FALSE <- function(use_gap) {
  if (is.null(use_gap) ||
      length(use_gap) != 1 ||
      !use_gap %in% c(TRUE, FALSE))
    stop("'is_clean' has to be either TRUE or FALSE")
}







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
                         .check_simple(any(floor(obj) != obj), argname, "has to be integer")
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
  .check_isnt_missing(obj)
  if (!allow_null) .check_isnt_null(obj, argname)
  else if (!is.null(obj)) {
    allowed <- c("ami", "nuc", if (allow_unt) "unt")
    .check_simple(obj %in% allowed, argname, paste0("has to be one of '", allowed, "'", collapse = "', '")) 
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
  .check_simple(length(obj_1) != length(obj_2), paste0(name_1, " and ", name_2), "have to have equal lenghts")
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
      stop("there are letters in alphabet that aren't suitable for given 'type' or 'is_clean' parameters")
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
    stop("all names of 'encoding' has to be letters from alphabet (elements of 'alphabet' attribute of 'sq')")
}
  
.check_has_both_UT <- function(has_U, has_T) {
  if (has_U && has_T)
    stop("'nucsq' sequences contains both 'U' and 'T' letters - should contain only one of them")
}

.check_has_A_no_UT <- function(has_U, has_T, has_A) {
  if (!has_U && !has_T && has_A)
    stop("'nucsq' sequences contains 'A' elements, but does not contain 'T' nor 'U' - unable to guess if it's dna or rna")
}

.check_is_dna_matches <- function(is_dna, has_U, has_T) {
  if ((is_dna && has_U) || (!is_dna && has_T)) 
    stop("if 'is_dna' is TRUE, sequences cannot contain 'U'; if is FALSE, sequences cannot contain 'T'")
}

.check_motifs_proper_alph <- function(motifs, type, alph = NULL) {
  if (type %in% c("ami", "nuc")) {
    if (!all(unlist(strsplit(motifs, "")) %in% c(.get_standard_alph(type, FALSE), "^", "$"))) 
      stop("motifs that you're searching for in the 'sq' object needs to consist of letters from its alphabet and optionally '^' or '$' characters")
  } else if (any(alph %in% c("^", "$", "?", "(", "=", ")", "\\", ".", "|", "+", "*", "{", "}", "[", "]"))) 
    stop("you cannot search for motifs if any of those characters: ^$?=()\\.|+*{}[] are elements of 'sq' alphabet; if you want to use them, please substitute those letters with some other using 'substitute_letters'")
}

.check_all_up_alph_proper <- function(up_alph, dest_alph) {
  if (!all(up_alph %in% dest_alph)) 
    stop("some sequences have levels that are invalid for given 'dest_type'; you can check them with 'get_invalid_letters' function and fix them with 'substitute_letters'")
}