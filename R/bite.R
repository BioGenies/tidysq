#' @export
bite <- function(x, index) {
  validate_sqtibble(x)
  if (!(is.numeric(index) && 
        floor(index) == index)) {
    stop("'index' has to be an integer vector")
  }
  x$sq <- lapply(x$sq, function(sq) sq[index])
  class(x$sq) <- "sqcol"
  x
}
