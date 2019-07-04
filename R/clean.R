#' @export
clean <- function(sqtbl, only_elements = FALSE) {
  #this function for sure can be vastly improved
  validate_sqtibble(sqtbl)
  types <- extract_sq_types(sqtbl)

  sqcol <- sqtbl[["sq"]]
  
  if (!only_elements) {
    inds_remove <- (sapply(sqcol, function(sq) any(sq %in% aminoacids_df[aminoacids_df[["amb"]], "one"])) & types == "aa") |
      (sapply(sqcol, function(sq) any(sq %in% nucleotides_df[nucleotides_df[["amb"]], "one"])) & types == "nuc") |
      (sapply(sqcol, function(sq) any(is.na(sq))) & types != "unt")
    sqtbl <- sqtbl[!inds_remove, ]
  } else {
    sqtbl[["sq"]][types == "aa"] <- lapply(sqcol[types == "aa"], function(sq) sq[!(sq %in% aminoacids_df[aminoacids_df[["amb"]], "one"])])
    sqtbl[["sq"]][types == "nuc"] <- lapply(sqcol[types == "nuc"], function(sq) sq[!(sq %in% nucleotides_df[nucleotides_df[["amb"]], "one"])])
    sqtbl[["sq"]][types != "unt"] <- lapply(sqcol[types != "unt"], function(sq) sq[!is.na(sq)])
  }
  if (any(types == "unt")) {
    warning("there are 'unt' sequences in sqtibble - they weren't changed")
  }
  sqtbl
}