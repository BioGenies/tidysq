#' @export
bite <- function(sqtbl, indices) {
  validate_sqtibble(sqtbl)

  if (!(is.numeric(indices) && 
        floor(indices) == indices)) {
    stop("'index' has to be an integer vector")
  }
  
  had_na <- any(sapply(sqtbl[["sq"]], function(sq) any(is.na(sq))))
  sqtbl[["sq"]] <- lapply(sqtbl[["sq"]], function(sq) sq[indices])
  has_na <- any(sapply(sqtbl[["sq"]], function(sq) any(is.na(sq))))
  if (has_na & !had_na) {
    warning("some sequences are subsetted with index bigger than lengh - NA's introduced")
  }
  set_sqcol(sqtbl)
}
