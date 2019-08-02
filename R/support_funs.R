#' @export
get_sq_alphabet <- function(sq) {
  validate_sq(sq)
  .get_alph(sq)
}

#' @export
get_sq_type <- function(sq) {
  validate_sq(sq)
  .get_sq_type(sq)
}

#' @export
get_sq_lengths <- function(sq) {
  validate_sq(sq)
  .get_lens(sq)
}

#' @export
is_null_sq <- function(sq) {
  validate_sq(sq)
  unlist(.apply_sq(sq, "char", "none", function(s) length(s) == 0 || 
                     length(s) == 1 && identical(s, "")))
}

#' @export
get_invalid_letters <- function(sq, dest_type) {
  validate_sq(sq)
  if (missing(dest_type) || 
     !(dest_type %in% c("ami", "nuc"))) {
    stop("'dest_type' should be either 'ami' or 'nuc'")
  }
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_char <- .get_na_char()
  dest_alph <- .get_standard_alph(dest_type, FALSE)
  dest_alph <- c(dest_alph, tolower(dest_alph))
  
  .apply_sq(sq, "char", "none", function(s) setdiff(s, dest_alph))
}