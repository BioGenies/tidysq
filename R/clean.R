#' @exportClass clnsq
#' @export
clean <- function(sq, only_elements = FALSE) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  is_clean <- .is_cleaned(sq)

  if (!is.logical(only_elements) ||
      is.na(only_elements) ||
      length(only_elements) != 1) {
    stop("'only_elements' has to be either TRUE or FALSE")
  }
  if (!(type %in% c("ami", "nuc"))) {
    stop("function 'clean' is meant to be used only with 'ami' or 'unt' sequences")
  }
  if (is_clean) {
    return(sq)
  }
  
  alph <- .get_alph(sq)
  alph_cln <- .get_standard_alph(type, TRUE)
  
  
  if (only_elements) {
    ret <- .apply_sq(sq, "char", "none", function(s) {
      pack_chars(s[(s %in% alph_cln)], alph_cln)
    })
  } else {
    ret <- .apply_sq(sq, "char", "none", function(s) {
      if (!all(s %in% alph_cln)) raw(1) else
        pack_chars(s, alph_cln)
    }) 
  }

  class(ret) <- c("clnsq", class(sq))
  .set_alph(ret, alph_cln)
}