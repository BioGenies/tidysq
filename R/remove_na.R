#' @export
remove_na <- function(sq, only_elements = FALSE) {
  validate_sq(sq)
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  
  if (only_elements) {
    ret <- lapply(sq, function(s) {
      s <- .bit_to_int(s, alph_size)
      .int_to_bit(s[s != na_val], alph_size)
    })
  } else {
    ret <- lapply(sq, function(s) {
      st <- .bit_to_int(s, alph_size)
      if (any(st == na_val)) raw(1) else s
    })
  }
  
  .set_class_alph(ret, sq)
}