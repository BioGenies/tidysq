#' @export
reverse <- function(sq) {
  validate_sq(sq)
  ret <- lapply(sq, rev)
  .set_class_alph(ret, sq)
}