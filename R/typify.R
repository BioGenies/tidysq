#' @export
typify <- function(sq, dest_type) {
  validate_sq(sq)
  if (missing(dest_type) || 
      !(dest_type %in% c("ami", "nuc"))) {
    stop("'dest_type' should be either 'ami' or 'nuc'")
  }
  
  type <- .get_sq_type(sq)
  if (type %in% c('ami', 'nuc')) {
    return(sq)
  }
  
  alph <- .get_alph(sq)
  up_alph <- unique(toupper(alph))
  dest_alph <- .get_standard_alph(dest_type, FALSE)
  
  if (!all(up_alph %in% dest_alph)) {
    stop("some sequences have levels that are invalid for given 'dest_type'; you can check them with 'get_invalid_letters' function and fix them with 'substitute_letters'")
  }
  
  if (!(length(alph) == length(up_alph))) {
    .handle_opt_txt("tidysq_typify_small_cap_let",
                    "in 'alphabet' attribute of 'sq' some letters show up as both lower and capital")
  }
  
  ret <- .apply_sq(sq, "char", "none", function(s) {
    s <- toupper(s)
    pack_chars(s, dest_alph)
  })
  
  ret <- .set_alph(ret, dest_alph)
  .set_class(ret, dest_type)
}