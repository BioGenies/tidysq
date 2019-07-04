#' @export
bite <- function(sqtbl, indices) {
  validate_sqtibble(sqtbl)
  if (!(is.numeric(indices) && 
        floor(indices) == indices)) {
    stop("'index' has to be an integer vector")
  }
  sqtbl[["sq"]] <- lapply(sqtbl[["sq"]], function(sq) sq[indices])
  set_sqcol(sqtbl)
}
