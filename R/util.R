#this functions are internal (at least by now) and are used only on sqtibbles we know are corect
#so theres no need to validate

set_sqcol <- function(sqtbl) {
  class(sqtbl[["sq"]]) <- "sqcol"
  sqtbl
}

extract_sq_types <- function(sqtbl) {
  sqclasses <- sapply(sqtbl[["sq"]], function(sq) intersect(class(sq), c("aasq", "nucsq", "untsq", "simsq")))
  dict <- c(aasq = "aa", nucsq = "nuc", untsq = "unt", simsq = "sim")
  dict[sqclasses]
}

set_clean <- function(sqtbl) {
  sqtbl[["sq"]] <- lapply(sqtbl[["sq"]], function(sq) {
    if (!("clnsq" %in% class(sq)) &&
        any(c("aasq", "nucsq") %in% class(sq))) {
      class(sq) <- c("clnsq", class(sq))
      sq
    } else {
      sq
    }
  })
  set_sqcol(sqtbl)
}

extract_is_clean <- function(sqtbl) {
  sapply(sqtbl[["sq"]], function(sq) "clnsq" %in% class(sq))
}