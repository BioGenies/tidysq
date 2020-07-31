#' @importFrom crayon green
#' @importFrom crayon silver
#' @export
format.sq <- function(x, ...) {
  max_sequences <- getOption("tidysq_p_max_sequences")
  use_color <- getOption("tidysq_p_use_color")
  
  rep("hidden", length(x))
}

#' @importFrom cli cli_text style_bold
#' @export
obj_print_header.sq <- function(x, ...) {
  if (length(.get_sq_type(x)) != 1)
    cli_text("sq (improper subtype!)")
  else
    cli_text("{vec_ptype_full(x)} sequences list")
}

#' @export
obj_print_data.sq <- function(x, ...) {
  # for some unknown reason had to copy this code from original obj_print_data
  if (length(x) == 0) {
    return()
  }
  print(format(x, ...), quote = FALSE)
  invisible(x)
}

#' @export
obj_print_footer.sq <- function(x, ..., num_lines = 0) {
  if (length(x) > num_lines)
    cat("printed ", num_lines, " out of ", length(x), sep = "")
  invisible(x)
}
