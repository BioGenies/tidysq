.validate_sq <- function(object, type = NULL) {
  argname <- deparse(substitute(object))
  if (!"sq" %in% class(object))
    stop(argname, " doesn't inherit class 'sq'", call. = FALSE)
  sqtype <- .get_sq_subclass(object)
  if (length(sqtype) != 1)
    stop("'object' should have exactly one of classes: 'amisq', 'nucsq', 'dnasq', 'rnasq', 'untsq', 'encsq', 'atpsq'", call. = FALSE)
  if (is.null(attr(object, "alphabet")))
    stop("'object' doesn't have 'alphabet' attribute", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.character(alph) &&
    !is.numeric(alph))
    stop("attribute 'alphabet' is neither a character nor a numeric vector", call. = FALSE)
  if (!is.list(object))
    stop("'object' isn't a list", call. = FALSE)
  if (!all(sapply(object, is.raw)))
    stop("'object' isn't a list of raw vectors", call. = FALSE)
  if (!is.null(type)) {
    switch(type,
           ami = .validate_amisq(object),
           nuc = .validate_nucsq(object),
           dna = .validate_dnasq(object),
           rna = .validate_rnasq(object),
           unt = .validate_untsq(object),
           atp = .validate_atpsq(object),
           enc = .validate_encsq(object),
           invisible(object)
    )
  }
}

.validate_dnasq <- function(object) {
  if (!"dnasq" %in% class(object))
    stop("'object' doesn't inherit class 'dnasq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!"clnsq" %in% class(object)) {
    if (!identical(alph, .get_standard_alph("dna", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard DNA alphabet", call. = FALSE)
  } else if (!identical(alph, .get_standard_alph("dna", TRUE)))
    stop("attribute 'alphabet' isn't identical to cleaned DNA alphabet", call. = FALSE)
  
  invisible(object)
}

.validate_rnasq <- function(object) {
  if (!"rnasq" %in% class(object))
    stop("'object' doesn't inherit class 'rnasq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!"clnsq" %in% class(object)) {
    if (!identical(alph, .get_standard_alph("rna", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard RNA alphabet", call. = FALSE)
  } else if (!identical(alph, .get_standard_alph("rna", TRUE)))
    stop("attribute 'alphabet' isn't identical to cleaned RNA alphabet", call. = FALSE)
  
  invisible(object)
}

.validate_nucsq <- function(object) {
  if (!"nucsq" %in% class(object))
    stop("'object' doesn't inherit class 'nucsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!"clnsq" %in% class(object)) {
    if (!identical(alph, .get_standard_alph("nuc", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard nucleotides alphabet", call. = FALSE)
  } else if (!identical(alph, .get_standard_alph("nuc", TRUE)))
    stop("attribute 'alphabet' isn't identical to cleaned nucleotides alphabet", call. = FALSE)
  
  invisible(object)
}

.validate_amisq <- function(object) {
  if (!"amisq" %in% class(object))
    stop("'object' doesn't inherit class 'amisq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!"clnsq" %in% class(object)) {
    if (!identical(alph, .get_standard_alph("ami", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard aminoacids alphabet", call. = FALSE)
  } else if (!identical(alph, .get_standard_alph("ami", TRUE)))
    stop("attribute 'alphabet' isn't identical to cleaned aminoacids alphabet", call. = FALSE)


  invisible(object)
}

.validate_untsq <- function(object) {
  if (!"untsq" %in% class(object))
    stop("'object' doesn't inherit class 'untsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector", call. = FALSE)

  invisible(object)
}

.validate_atpsq <- function(object) {
  if (!"atpsq" %in% class(object))
    stop("'object' doesn't inherit class 'atpsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector", call. = FALSE)

  invisible(object)
}

.validate_encsq <- function(object) {
  if (!"encsq" %in% class(object))
    stop("'object' doesn't inherit class 'encsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.numeric(alph))
    stop("attribute 'alphabet' isn't a numeric vector", call. = FALSE)

  invisible(object)
}
