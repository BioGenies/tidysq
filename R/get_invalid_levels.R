#' @export
get_invalid_levels <- function(sqtbl, dest_type, only_levels = TRUE) {
  validate_sqtibble(sqtbl)
  if (missing(dest_type) || 
     !(dest_type %in% c("aa", "nuc"))) {
    stop("'dest_type' should be either 'aa' or 'nuc'")
  }
  if (!is.logical(only_levels) ||
      (length(only_levels) != 1) ||
      is.na(only_levels)) {
    stop("'only_levels' has to be either TRUE or FALSE")
  }
  
  sqtypes <- extract_sq_types(sqtbl)
  
  if (!all(sqtypes %in% c("unt", dest_type))) {
    stop("each sequence in 'sqtbl' should have either 'unt' type or type given in 'dest_type'")
  }
  
  inv_levels <- extract_inv_lvls(sqtbl, dest_type, sqtypes)
  if (only_levels) {
    unique(unlist(inv_levels))
  } else {
    inv_levels
  }
}
