#this functions are internal (at least by now) and are used only on sqtibbles we know are corect
#so theres no need to validate

set_sqcol <- function(sqtbl) {
  class(sqtbl[["sq"]]) <- "sqcol"
  sqtbl
}

extract_sq_types <- function(sqtbl) {
  classes_list <- lapply(sqtbl[["sq"]], function(sq) intersect(class(sq), c("aasq", "nucsq", "untsq")))
  if (!all(sapply(classes_list, function(class_vec) length(class_vec) == 1))) {
    stop("not every sequence in sqtibble has appropriate type (each should have exactly one of subclasses: 'aasq', 'nucsq', 'untsq')")
  }
  classes_list <- unlist(classes_list)
  if (!all(classes_list %in% c("aasq", "nucsq", "untsq"))) {
    stop("not every sequence in sqtibble has appropriate type (each should have exactly one of subclasses: 'aasq', 'nucsq', 'untsq')")
  }
  dict <- c(aasq = "aa", nucsq = "nuc", untsq = "unt")
  dict[classes_list]
}