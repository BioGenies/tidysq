#' @export
remove_na <- function(sqtbl, only_elements = FALSE) {
  validate_sqtibble(sqtbl)
  
  sqcol <- sqtbl[["sq"]]
  if (!only_elements) {
    inds_remove <- sapply(sqcol, function(sq) any(is.na(sq)))
    sqtbl <- sqtbl[!inds_remove, ]
  } else {
    sqtbl[["sq"]] <- lapply(sqcol, function(sq) sq[!is.na(sq)])
    sqtbl <- set_sqcol(sqtbl)
  }
  sqtbl
}