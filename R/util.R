#this functions are internal (at least by now) and are used only on sq objects that we know are corect
#so there's no need to validate

.get_alph_size <- function(alph) {
  ceiling(log2(length(alph) + 2))
}

.get_na_val <- function(alph) {
  2 ^ .get_alph_size(alph) - 1
}

.get_alph <- function(sq) {
  attr(sq, "alphabet")
}

.get_real_alph <- function(sq) {
  unique(unlist(strsplit(sq, "")))
}

.get_sq_subclass <- function(sq) {
  intersect(class(sq), c("amisq", "nucsq", "untsq", "simsq", "atpsq", "encsq"))
}

.get_sq_type <- function(sq) {
  sqclasses <- intersect(class(sq), c("amisq", "nucsq", "untsq", "simsq", "atpsq", "encsq"))
  dict <- c(amisq = "ami", nucsq = "nuc", untsq = "unt", simsq = "sim", atpsq = "atp", encsq = "encsq")
  dict[sqclasses]
}

.get_standard_alph <- function(type, is_clean) {
       if (type == "ami" &&  is_clean) aminoacids_df[!aminoacids_df[["amb"]], "one"]
  else if (type == "ami" && !is_clean) aminoacids_df[, "one"]
  else if (type == "nuc" &&  is_clean) nucleotides_df[!nucleotides_df[["amb"]], "one"]
  else if (type == "nuc" && !is_clean) nucleotides_df[, "one"]
}

.is_cleaned <- function(sq) {
  "clnsq" %in% class(sq)
}

.set_class <- function(sq, type, is_clean = FALSE) {
  class(sq) <- c(if (is_clean) "clnsq" else NULL, paste0(type, "sq"), "sq")
  sq
}

.set_alph <- function(sq, alph) {
  attr(sq, "alphabet") <- alph
  sq
}


.set_class_alph <- function(new_sq, sq) {
  class(new_sq) <-class(sq)
  attr(new_sq, "alphabet") <- .get_alph(sq)
  new_sq
}

.construct_sq_s <- function(sq, alph, classes) {
  sq <- .bitify_sq(sq, alph)
  attr(sq, "alphabet") <- alph
  class(sq) <- classes
  sq
}

.guess_ami_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("ami", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("ami", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: aminoacids_df)")
}

.guess_nuc_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("nuc", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("nuc", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: nucleotides_df)")
}

.guess_sq_type <- function(sq) {
  real_alph <- toupper(.get_real_alph(sq))
  if (all(real_alph %in% .get_standard_alph("nuc", FALSE))) "nuc"
  else if (all(real_alph %in% .get_standard_alph("ami", FALSE))) "ami"
  else "unt"
}