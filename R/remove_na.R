#' @export
remove_na <- function(sq, only_elements = FALSE) {
  validate_sq(sq)
  
  if (only_elements) {
    ret <- lapply(sq, function(s) s[!is.na(s)])
  } else {
    ret <- lapply(sq, function(s) if (any(is.na(s))) integer(0) else s)
  }
  
  .set_class_alph(ret, sq)
}

#' @exportMethod na.omit sq
#' @export
na.omit.sq <- function(object, ...) {
  remove_na(object, ...)
}