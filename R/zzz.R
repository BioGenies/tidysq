.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    tidysq_max_sq_print_width = 15
  )
  
  toset <- !(names(new_options) %in% names(prev_options))
  if (any(toset)) {
    options(new_options)
  }
  
  invisible()
}