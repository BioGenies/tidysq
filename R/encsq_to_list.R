#' @export
encsq_to_list <- function(encsq) {
  validate_sq(encsq, type = "enc")
  
  alph <- .get_alph(encsq)
  .apply_sq(encsq, "int", "none", function(s) alph[s])
}