.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    tidysq_a_bite_na = "warning",
    tidysq_a_cln_sub_letters = "warning",
    tidysq_a_no_given_enc = "warning",
    tidysq_a_typify_small_cap_let = "warning",
    tidysq_g_fast_mode = FALSE,
    tidysq_p_max_pillar_width = 15,
    tidysq_p_max_sequences = 10,
    tidysq_p_na_letter = "!",
    tidysq_p_use_color = TRUE,
    
    tidysq_NA_letter = "!",
    tidysq_safe_mode = FALSE
  )
  
  unset_inds <- !(names(new_options) %in% names(prev_options))
  if (any(unset_inds)) {
    options(new_options[unset_inds])
  }
  
  invisible()
}