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
    handle_opt_txt("tidysq_bite_na_action",
                   "some sequences are subsetted with index bigger than length - NA's introduced")
  }
  set_sqcol(sqtbl)
}
