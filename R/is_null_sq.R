#' @export

is_null_sq <- function(sq) {
  validate_sq(sq)
  lengths(sq) == 0
}