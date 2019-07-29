#' @export
get_invalid_letters <- function(sq, dest_type) {
  validate_sq(sq)
  if (missing(dest_type) || 
     !(dest_type %in% c("ami", "nuc"))) {
    stop("'dest_type' should be either 'ami' or 'nuc'")
  }
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  alph <- .rawize_alph(alph)
  na_char <- .get_na_char()
  dest_alph <- if (dest_type == "ami") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  dest_alph <- c(dest_alph, tolower(dest_alph))
  
  .apply_sq(sq, "chars", function(s) setdiff(s, dest_alph))
}
