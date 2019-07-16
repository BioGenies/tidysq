#' @export

is_null_sq <- function(sq) {
  validate_sq(sq)
  sapply(sq, function(s) s == raw(1))
}