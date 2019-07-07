#' @export
set_sq_types <- function(sqtbl, dest_type) {
  validate_sqtibble(sqtbl)
  if (missing(dest_type) || 
      !(dest_type %in% c("aa", "nuc"))) {
    stop("'dest_type' should be either 'aa' or 'nuc'")
  }
  
  sqtypes <- extract_sq_types(sqtbl)
  if (!all(sqtypes %in% c("unt", dest_type))) {
    stop("each sequence in 'sqtbl' should have either 'unt' type or type given in 'dest_type'")
  }
  
  inv_lvls_alph <- unique(unlist(extract_inv_lvls(sqtbl, dest_type, sqtypes)))
  if (length(inv_lvls_alph) > 0) {
    stop("some sequences have levels that are invalid for given 'dest_type'; you can check them with 'get_invalid_levels' function and fix them with 'substitute_invalid_levels' or 'drop_invalid_levels'")
  }
  
  dest_alph <- if (dest_type == "aa") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  sqtbl[["sq"]][sqtypes != dest_type] <- lapply(sqtbl[["sq"]][sqtypes != dest_type] , function(sq) {
    sq <- factor(as.character(sq), levels = dest_alph)
    class(sq) <- c(paste0(dest_type, "sq"), "sq", "factor")
    sq
  })
  sqtbl
}