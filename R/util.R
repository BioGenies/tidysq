#this functions are internal (at least by now) and are used only on sqtibbles we know are corect
#so theres no need to validate

set_sqcol <- function(sqtbl) {
  class(sqtbl[["sq"]]) <- "sqcol"
  sqtbl
}

extract_sq_types <- function(sqtbl) {
  classes_list <- sapply(sqtbl[["sq"]], function(sq) intersect(class(sq), c("aasq", "nucsq", "untsq", "simsq")))
  dict <- c(aasq = "aa", nucsq = "nuc", untsq = "unt", simsq = "sim")
  dict[classes_list]
}