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

.is_cleaned <- function(sq) {
  "clnsq" %in% class(sq)
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
  if (all(real_alph %in% aminoacids_df[!aminoacids_df[["amb"]], "one"]))
    TRUE
  else if (all(real_alph %in% aminoacids_df[, "one"]))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: aminoacids_df)")
}

.guess_nuc_is_clean <- function(real_alph) {
  if (all(real_alph %in% nucleotides_df[!nucleotides_df[["amb"]], "one"]))
    TRUE
  else if (all(real_alph %in% nucleotides_df[, "one"]))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: nucleotides_df)")
}

.guess_sq_type <- function(sq) {
  real_alph <- toupper(.get_real_alph(sq))
  if (all(real_alph %in% nucleotides_df[, "one"])) "nuc"
  else if (all(real_alph %in% aminoacids_df[, "one"])) "ami"
  else "unt"
}