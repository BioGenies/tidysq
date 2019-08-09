.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    tidysq_bite_na_action = "warning",
    tidysq_subsitute_letters_cln = "warning",
    tidysq_typify_small_cap_let = "warning",
    tidysq_encode_no_given_action = "warning",
    tidysq_max_pillar_sq_width = 15,
    tidysq_max_print_sequences = 10,
    tidysq_colorful_sq_print = TRUE,
    tidysq_na_print_char = "!",
    tidysq_no_check_mode = FALSE
  )
  
  unset_inds <- !(names(new_options) %in% names(prev_options))
  if (any(unset_inds)) {
    options(new_options[unset_inds])
  }
  
  invisible()
}