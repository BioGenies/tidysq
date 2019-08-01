#' @exportClass simsq
#' @export
simplify <- function(sq, encoding) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  
  if (!is.character(encoding) ||
      any(!(encoding %in% c(letters, "-")))) {
    stop("'encoding' should be a named vector, where names are upper latin letters or '-' and elements are lower latin letters and '-'")
  }
  
  alph <- .get_alph(sq)
  if (!all(alph %in% names(encoding))) {
    stop("each of letters (elements of 'alphabet' attribute of 'sq') should appear in names of encoding")
  }
  
  new_alph <- sort(unique(encoding))
  inds_fun <- match(encoding, new_alph)
  names(inds_fun) <- as.character(match(names(encoding), alph))
  
  new_alph_size <- .get_alph_size(new_alph)
  
  ret <- .apply_sq(sq, "int", "none", function(s) {
    pack_ints(inds_fun[s], new_alph_size)
  })

  ret <- .set_alph(ret, new_alph)
  .set_class(ret, "sim")
}

