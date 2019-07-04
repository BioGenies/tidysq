set_sqcol <- function(sqtbl) {
  class(sqtbl[["sq"]]) <- "sqcol"
  sqtbl
}

extract_sq_types <- function(sqtbl) {
  ret <- lapply(x[["sq"]], function(sq) intersect(class(sq), c("aa", "nuc", "unt")))
  if (!all(ret %in% c("aasq", "nucsq", "untsq"))) {
    stop("not every sequence in sqtibble has appropriate type (each should have exactly one of subclasses: 'aasq', 'nucsq', 'untsq')")
  }
  dict <- c(aasq = "aa", nucsq = "nuc", ambsq = "unt")
  dict[ret]
}