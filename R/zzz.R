.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    tidysq_bite_na_action = "warning",
    tidysq_subsitute_letters_cln = "warning",
    tidysq_typify_small_cap_let = "warning",
    tidysq_encode_nogap_action = "none",
    tidysq_encode_noname_action = "warning",
    tidysq_max_sq_print_width = 15,
    tidysq_colorful_sq_print = TRUE,
    tidysq_na_print_char = "!"
  )
  
  unset_inds <- !(names(new_options) %in% names(prev_options))
  if (any(unset_inds)) {
    options(new_options[unset_inds])
  }
  
  invisible()
}