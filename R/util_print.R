#' @importFrom cli cli_text
#' @export
obj_print_header.sq <- function(x, ...) {
  if (length(.get_sq_type(x)) != 1)
    cli_text("sq (improper subtype!)")
  else
    cli_text("{vec_ptype_full(x)} sequences list",
             if (length(x) == 0) " of length 0" else ":")
}

#' @importFrom cli cli_text
#' @export
obj_print_data.sq <- function(x, ...) {
  # for some unknown reason had to copy this code from original obj_print_data
  if (length(x) == 0) {
    return()
  }
  cli_text(format(x, ...))
  invisible(x)
}

#' @importFrom cli cli_text style_italic
#' @export
obj_print_footer.sq <- function(x, ...,
                                max_sequences = getOption("tidysq_p_max_sequences")) {
  if (length(x) > max_sequences)
    cli_text(style_italic("printed {max_sequences} out of {length(x)}"))
  invisible(x)
}
