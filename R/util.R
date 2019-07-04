set_sqcol <- function(sqtbl) {
  class(sqtbl[["sq"]]) <- "sqcol"
  sqtbl
}