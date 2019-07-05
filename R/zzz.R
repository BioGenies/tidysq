.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    tidysq_constr_mtype_action = "warning", 
    tidysq_constr_unt_action = "warning",
    tidysq_bite_na_action = "warning",
    tidysq_max_sq_print_width = 15,
    tidysq_colorful_sq_print = TRUE
  )
  
  unset_inds <- !(names(new_options) %in% names(prev_options))
  if (any(unset_inds)) {
    options(new_options[unset_inds])
  }
  
  invisible()
}