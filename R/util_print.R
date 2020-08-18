#' @importFrom cli cli_text
#' @export
obj_print_header.sq <- function(x, ...) {
  if (length(.get_sq_type(x)) != 1)
    cli_text("sq (improper subtype!)")
  else
    cli_text("{vec_ptype_full(x)} sequences list",
             if (length(x) == 0) " of length 0" else ":")
}

#' @importFrom cli cat_line
#' @export
obj_print_data.sq <- function(x, ...) {
  # for some unknown reason had to copy this code from original obj_print_data
  if (length(x) == 0) {
    return()
  }
  cat_line(format(x, ...))
  invisible(x)
}

#' @importFrom cli cli_text
#' @export
obj_print_footer.sq <- function(x, ...,
                                max_sequences = getOption("tidysq_p_max_sequences")) {
  if (length(x) > max_sequences)
    cli_text("printed {max_sequences} out of {length(x)}")
  invisible(x)
}

#' @importFrom cli col_green col_silver
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.sq <- function(x, ...) {
  # color NA's
  na_character(alphabet(x)) <- col_silver(.get_na_char())
  
  .pillar_shaft_sq(x, "", col_green)
}

#' @importFrom cli col_cyan
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.encsq <- function(x, ...) {
  alphabet(x) <- format(alphabet(x), digits = 1, scientific = FALSE)
  
  .pillar_shaft_sq(x, "\u00a0", col_cyan)
}
