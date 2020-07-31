#' @importFrom crayon green
#' @importFrom crayon silver
#' @export
format.sq <- function(x, ...) {
  max_sequences <- getOption("tidysq_p_max_sequences")
  use_color <- getOption("tidysq_p_use_color")
  
  rep("hidden", length(x))
}

#' @export
obj_print_data.sq <- function(x, ...) {
  # for some unknown reason had to copy this code from original obj_print_data
  if (length(x) == 0) {
    return()
  }
  print(format(x), quote = FALSE)
  invisible(x)
}

#' @export
obj_print_footer.sq <- function(x, ...) {
  cat("peek-a-boo!")
  invisible(x)
}
