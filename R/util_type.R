guess_sq_type <- function(x) {
  
}

interpret_type <- function(name) {
  
}

type_as_class <- function(type) {
  paste0("sq_", type)
}

# TODO: sort that rubbish

.guess_sq_type_subtype <- function(sq) {
  real_alph <- toupper(.get_real_alph(sq))
  .guess_type_subtype_by_alph(real_alph)
}

.guess_type_subtype_by_alph <- function(alph) {
  # TODO: any better idea for it?
  possib_rets <- list(list(type = "dna", is_clean = TRUE),
                      list(type = "rna", is_clean = TRUE),
                      list(type = "ami", is_clean = TRUE),
                      list(type = "dna", is_clean = FALSE),
                      list(type = "rna", is_clean = FALSE),
                      list(type = "ami", is_clean = FALSE))
  for (ret in possib_rets)
    if (all(alph %in% .get_standard_alph(ret[["type"]], ret[["is_clean"]]))) return(ret)
  list(type = "unt", is_clean = NULL)  
}

.get_sq_type <- function(sq) {
  sqclasses <- intersect(class(sq), c("amisq", "dnasq", "rnasq", "untsq", "atpsq", "encsq"))
  dict <- c(amisq = "ami", dnasq = "dna", rnasq = "rna", untsq = "unt", atpsq = "atp", encsq = "enc")
  dict[sqclasses]
}
