.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    tidysq_NA_letter = "!",
    tidysq_safe_mode = FALSE,
    tidysq_on_warning = "warning",
    tidysq_pillar_max_width = 15,
    tidysq_print_max_sequences = 10,
    tidysq_print_use_color = TRUE
  )
  
  unset_inds <- !(names(new_options) %in% names(prev_options))
  if (any(unset_inds)) {
    options(new_options[unset_inds])
  }
  
  invisible()
}