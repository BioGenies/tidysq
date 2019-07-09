#this functions are internal (at least by now) and are used only on sqtibbles we know are corect
#so theres no need to validate

.get_alph <- function(sq) {
  attr(sq, "alphabet")
}

.get_sq_subclass <- function(sq) {
  intersect(class(sq), c("amisq", "nucsq", "untsq", "simsq", "atpsq"))
}

.get_sq_type <- function(sq) {
  sqclasses <- intersect(class(sq), c("amisq", "nucsq", "untsq", "simsq", "atpsq"))
  dict <- c(amisq = "ami", nucsq = "nuc", untsq = "unt", simsq = "sim", atpsq = "atp")
  dict[sqclasses]
}

.is_cleaned <- function(sq) {
  "clnsq" %in% class(sq)
}

#####

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

extract_inv_lvls <- function(sqtbl, dest_type, sqtypes) {
  inv_levels <- vector("list", nrow(sqtbl))
  dest_alph <- if (dest_type == "aa") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  
  inv_levels[sqtypes != dest_type] <- lapply(sqtbl[["sq"]][sqtypes != dest_type], 
                                             function(sq) setdiff(levels(sq), dest_alph))
}