#' @export

is_null_sq <- function(sq) {
  validate_sq(sq)
  unlist(.apply_sq(sq, "char", "none", function(s) length(s) == 0 || 
                     length(s) == 1 && identical(s, "")))
}