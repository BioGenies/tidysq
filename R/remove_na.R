#' @export
remove_na <- function(sq, only_elements = FALSE) {
  validate_sq(sq)
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  
  if (only_elements) {
    ret <- .apply_sq(sq, "int", "int", function(s) {
      s[s != na_val]
    }) 
  } else {
    ret <- lapply(sq, function(s) {
      st <- unpack_ints(s, alph_size)
      if (any(st == na_val)) raw(1) else s
    })
  }
  
  .set_class_alph(ret, sq)
}