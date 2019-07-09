#' @export
get_invalid_levels <- function(sq, dest_type) {
  validate_sq(sq)
  if (missing(dest_type) || 
     !(dest_type %in% c("ami", "nuc"))) {
    stop("'dest_type' should be either 'ami' or 'nuc'")
  }
  
  alph <- .get_alph(sq)
  dest_alph <- if (dest_type == "ami") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  dest_alph <- c(dest_alph, tolower(dest_alph))
  
  inv_levels <- lapply(sq, function(s) setdiff(alph[s], dest_alph))
  inv_levels
}
