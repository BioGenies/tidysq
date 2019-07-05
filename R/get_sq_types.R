#' @export
get_sq_types <- function(sqtbl) {
  validate_tibble(sqtbl)
  ret <- extract_sq_types(sqtbl)
  ret <- enframe(ret, name = NULL)
  names(ret) <- "type"
  ret
}