#' @export
clean <- function(x, only_elements = FALSE) {
  types <- extract_sq_types(x)
  
  if (!only_elements) {
    indicies_remove <- (sapply(unclass(x[["sq"]]), function(sq) any(sq %in% aminoacids_df[aminoacids_df[["amb"]], "one"])) & types == "aa") |
                       (sapply(unclass(x[["sq"]]), function(sq) any(sq %in% nucleotides_df[nucleotides_df[["amb"]], "one"])) & types == "nuc") |
                       (sapply(unclass(x[["sq"]]), function(sq) any(is.na(sq)) & types != "unt"))
    if (any(types == "unt")) {
      warning("there are 'unt' sequences in sqtibble - they weren't removed")
    }
    x[!indicies_remove,]
  } else {
    x[["sq"]][types == "aa"] <- lapply(x[["sq"]][types == "aa"], function(sq) sq[!(sq %in% aminoacids_df[aminoacids_df[["amb"]], "one"])])
    x[["sq"]][types == "nuc"] <- lapply(x[["sq"]][types == "nuc"], function(sq) sq[!(sq %in% nucleotides_df[nucleotides_df[["amb"]], "one"])])
    x[["sq"]][types != "unt"] <- lapply(x[["sq"]][types != "unt"], function(sq) sq[!is.na(sq)])
    if (any(types == "unt")) {
      warning("there are 'unt' sequences in sqtibble - they weren't changed")
    }
    x
  }
}

extract_sq_types <- function(x) {
  validate_sqtibble(x)
  ret <- sapply(x[["sq"]], function(sq) class(sq)[1])
  if (!all(ret %in% c("aasq", "nucsq", "untsq"))) {
    stop("not every sequence has a type (each should have subclass 'aasq', 'nucsq', 'untsq')")
  }
  dict <- c(aasq = "aa", nucsq = "nuc", ambsq = "unt")
  dict[ret]
}
