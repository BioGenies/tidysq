#' @export
drop_invalid_levels <- function(sqtbl, dest_type, replacement = NA) {
  validate_sqtibble(sqtbl)
  if (missing(dest_type) || 
      !(dest_type %in% c("aa", "nuc"))) {
    stop("'dest_type' should be either 'aa' or 'nuc'")
  }
  
  any_char <- if (dest_type == "aa") "X" else "N"
  if (!(replacement %in% c(NA, any_char))) {
    stop("'replacement' has to be either NA or 'X' when 'dest_type' is 'aa' or 'N' when 'dest_type' is 'nuc'")
  }
  
  sqtypes <- extract_sq_types(sqtbl)
  if (!all(sqtypes %in% c("unt", dest_type))) {
    stop("each sequence in 'sqtbl' should have either 'unt' type or type given in 'dest_type'")
  }
  
  inv_lvls_alph <- unique(unlist(extract_inv_lvls(sqtbl, dest_type, sqtypes)))
  dest_alph <- if (dest_type == "aa") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  
  sqtbl[["sq"]][sqtypes != dest_type] <- lapply(sqtbl[["sq"]][sqtypes != dest_type] , function(sq) {
    lvls <- levels(sq)
    lvls[!(lvls %in% dest_alph)] <- replacement
    levels(sq) <- lvls 
    sq
  })
  sqtbl
  
}