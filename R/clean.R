#' @exportClass clnsq
#' @export
clean <- function(sqtbl, only_elements = FALSE) {
  #this function for sure can be written way more efficiently
  validate_sqtibble(sqtbl)
  types <- extract_sq_types(sqtbl)
  ind_clean <- extract_is_clean(sqtbl)

  sqcol <- sqtbl[["sq"]]
  
  if (!only_elements) {
    inds_remove <- (sapply(sqcol, function(sq) any(sq %in% aminoacids_df[aminoacids_df[["amb"]], "one"])) & types == "aa") |
      (sapply(sqcol, function(sq) any(sq %in% nucleotides_df[nucleotides_df[["amb"]], "one"])) & types == "nuc")
    sqtbl <- sqtbl[!inds_remove, ]
    types <- types[!inds_remove]
    sqtbl[["sq"]][types == "aa"] <- lapply(sqtbl[["sq"]][types == "aa"], function(sq) {
      sq <- factor(as.character(sq), levels = aminoacids_df[!aminoacids_df[["amb"]], "one"])
      class(sq) <- c("aasq", "sq", "factor")
      sq
    })
    sqtbl[["sq"]][types == "nuc"] <- lapply(sqtbl[["sq"]][types == "nuc"], function(sq) {
      sq <- factor(as.character(sq), levels = nucleotides_df[!nucleotides_df[["amb"]], "one"])
      class(sq) <- c("nucsq", "sq", "factor")
      sq
    })
    
  } else {
    sqtbl[["sq"]][types == "aa" & !ind_clean] <- lapply(
      sqcol[types == "aa"], function(sq) {
        sq <- sq[!(sq %in% aminoacids_df[aminoacids_df[["amb"]], "one"])]
        sq <- factor(as.character(sq), levels = aminoacids_df[!aminoacids_df[["amb"]], "one"])
        class(sq) <- c("aasq", "sq", "factor")
        sq
      })
    sqtbl[["sq"]][types == "nuc" & !ind_clean] <- lapply(
      sqcol[types == "nuc"], function(sq) {
        sq <- sq[!(sq %in% nucleotides_df[nucleotides_df[["amb"]], "one"])]
        sq <- factor(as.character(sq), levels = nucleotides_df[!nucleotides_df[["amb"]], "one"])
        class(sq) <- c("nucsq", "sq", "factor")
        sq
      })
  }
  
  if (any(types == "unt")) {
    warning("there are 'unt' or 'sim' sequences in sqtibble - they weren't changed")
  }
  set_clean(sqtbl)
}