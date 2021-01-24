#' @export
paste <- function(...)
  UseMethod("paste")

#' @export
paste.default <- base::paste

#' @family order_functions
#' @export
paste.sq <- function(...,
                     NA_letter = getOption("tidysq_NA_letter")) {
  # Throws error when there is no common size
  vec_size_common(...)
  assert_string(NA_letter, min.chars = 1)

  CPP_paste(vec_cast_common(...), NA_letter)
}
