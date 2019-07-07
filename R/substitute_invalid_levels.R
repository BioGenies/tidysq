#' @export
substitue_invalid_levels <- function(sqtbl, dest_type, encoding) {
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
  if (!all(names(encoding) %in% inv_lvls_alph)) {
    stop("all names of 'encoding' has to be levels that are invalid")
  }
  dest_alph <- if (dest_type == "aa") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  encoding <- toupper(encoding)
  if (!all(encoding %in% dest_alph)) {
    stop("all elements of 'encoding' should be elements of destination alphabet")
  }
  
  sqtbl[["sq"]][sqtypes != dest_type] <- lapply(sqtbl[["sq"]][sqtypes != dest_type] , function(sq) {
    lvls <- levels(sq)
    matched <- match(names(encoding), lvls)
    lvls[matched[!is.na(matched)]] <- encoding[names(encoding)[!is.na(matched)]]
    levels(sq) <- lvls 
    sq
  })
  sqtbl
}