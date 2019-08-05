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

.check_type_or_nonst_alph <- function(type, is_clean, non_standard) {
  if (!(is.null(type) && 
        is.null(is_clean)) &&
      !is.null(non_standard))
    stop("if you specify 'non_standard', you cannot specify neither 'type' nor 'is_clean'")
}

.check_nonst_proper_char <- function(non_standard) {
  if (!is.character(sq) ||
      length(sq) == 0 ||
      any(is.na(sq)) ||
      any(is.null(sq))) 
    stop("'non_standard' has to be a positive-length character vector without NA")
}

.check_nonst_nchar <- function(non_standard) {
  if (any(nchar(non_standard) < 2)) 
    stop("non standard letters specified in 'non_standard' parameter have to have more than one character")
}

.check_alph_matches_type <- function(alph, type, is_clean) {
  if (!is.null(type) && !(type == "unt")) {
    alph <- toupper(alph)
    if (is.null(is_clean)) is_clean <- FALSE
    if (!all(alph %in% .get_standard_alph(type, is_clean))) {
      stop("there are letters in alphabet in file that aren't suit for given 'type' or 'is_clean' parameters")
    }
  }
}

.check_alph_length <- function(alph) {
  if (length(alph) > 30) 
    stop("max length of alphabet is 30 letters")
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
      is.nan(nchar)|| 
      !is.finite(nchar) ||
      nchar <= 0) 
    stop("'nchar' has to be positive integer indicating max number of elements of sequence in single line in file")
}

.check_eq_lens <- function(sq, name) {
  if (length(sq) != length(name)) 
    sstop("'name' has to have length equal to 'sq'")
}