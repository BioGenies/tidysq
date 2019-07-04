#' @export
reverse <- function(sqtbl) {
  validate_sqtibble(sqtbl)
  sqtbl[["sq"]] <- lapply(sqtbl[["sq"]], rev)
  set_sqcol(sqtbl)
}