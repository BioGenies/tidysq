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