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

.set_class_alph <- function(new_sq, sq) {
  class(new_sq) <-class(sq)
  attr(new_sq, "alphabet") <- .get_alph(sq)
  new_sq
}