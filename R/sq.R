#' @exportClass sq
#' @export
construct_sq <- function(sq, type = "unt") {
  if (!type %in% c("unt", "ami", "nuc")) {
    stop("'type' has to be one of 'unt', 'ami', 'nuc'")
  }
  if (!is.character(sq)) {
    stop("'sq' has to be a vector of strings", call. = FALSE)
  }
  
  if (type == "ami") {
    construct_amisq(sq)
  } else if (type =="nuc") {
    construct_nucsq(sq)
  } else {
    construct_untsq(sq)
  }
}

#' @exportClass nucsq
construct_nucsq <- function(sq) {
  sq <- toupper(sq)
  alph <- nucleotides_df[,"one"]
  sq <- strsplit(sq, "")
  is_nuc_sq <- all(unlist(sq) %in% alph)
  if (!is_nuc_sq) {
    stop("each of letters should be in nucleotide alphabet (one of nucleotides_df[,'one'])")
  }
  
  object <- lapply(sq, function(s) match(s, alph))
  attr(object, "alphabet") <- alph
  class(object) <- c("nucsq", "sq")
  object
}

#' @exportClass amisq
construct_amisq <- function(sq) {
  sq <- toupper(sq)
  alph <- aminoacids_df[,"one"]
  sq <- strsplit(sq, "")
  is_ami_sq <- all(unlist(sq) %in% alph)
  if (!is_ami_sq) {
    stop("each of letters should be in aminoacids alphabet (one of aminoacids_df[,'one'])")
  }
  
  object <- lapply(sq, function(s) match(s, alph))
  attr(object, "alphabet") <- alph
  class(object) <- c("amisq", "sq")
  object
}

#' @exportClass untsq
construct_untsq <- function(sq) {
  sq <- strsplit(sq, "")
  alph <- unique(unlist(sq))

  object <- lapply(sq, function(s) match(s, alph))
  attr(object, "alphabet") <- alph
  class(object) <- c("untsq", "sq")
  object
}

#'
validate_sq <- function(object, type = NULL) {
  if (!"sq" %in% class(object)) {
    stop("'object' doesn't inherit class 'sq'")
  } 
  sqtype <- .get_sq_subclass(object)
  if (!(length(sqtype) == 1)) {
    stop("'object' should have exactly one of types: 'ami', 'nuc', 'unt', 'sim', 'atp")
  }
  if (is.null(attr(object, "alphabet"))) {
    stop("'object' doesn't 'alphabet' attribute")
  }
  if (is.null(attr(object, "alphabet"))) {
    stop("'object' doesn't 'alphabet' attribute")
  }
  alph <- .get_alph(object)
  if (!is.character(alph)) {
    stop("attribute 'alphabet' isn't character vector")
  }
  #assumption about length of one of each character - this can be changed in future
  if (!all(sapply(alph, length) == 1)) {
    stop("attribute 'alphabet' have elements, that aren't one character long")
  }
  if (!is.list(object)) {
    stop("'object' isn't a list")
  }
  if (!all(sapply(object, is.integer))) {
    stop("'object' isn't a list of integer vectors")
  }
  if (max(unlist(object), na.rm = TRUE) > length(alph)) {
    stop("'alphabet' attribute has less elements than are different values in 'object'")
  }
  if (!is.null(type)) {
    switch (type,
      ami = validate_amisq(object),
      nuc = validate_nucsq(object),
      unt = validate_untsq(object),
      sim = validate_simsq(object),
      atp = validate_atpsq(object)
    )
  }
  invisible(object)
}

#'
validate_nucsq <- function(object) {
  validate_sq(object)
  if (!"nucsq" %in% class(object)) {
    stop("'object' doesn't inherit class 'nucsq'")
  } 
  alph <- .get_alph(object)
  if (!("clnsq" %in% class(object))) {
    if (!identical(alph, nucleotides_df[,"one"])) {
      stop("attribute 'alphabet' isn't identical to standard nucleotides alphabet")
    }
  } else {
    if (!identical(alph, nucleotides_df[!nucleotides_df[["amb"]],"one"])) {
      stop("attribute 'alphabet' isn't identical to cleaned nucleotides alphabet")
    }
  }
  invisible(object)
}

#'
validate_amisq <- function(object) {
  validate_sq(object)
  if (!"amisq" %in% class(object)) {
    stop("'object' doesn't inherit class 'amisq'")
  } 
  alph <- .get_alph(object)
  if (!("clnsq" %in% class(object))) {
    if (!identical(alph, aminoacids_df[,"one"])) {
      stop("attribute 'alphabet' isn't identical to standard aminoacids alphabet")
    }
  } else {
    if (!identical(alph, aminoacids_df[!aminoacids_df[["amb"]],"one"])) {
      stop("attribute 'alphabet' isn't identical to cleaned aminoacids alphabet")
    }
  }
  
  invisible(object)
}

#'
validate_untsq <- function(object) {
  validate_sq(object)
  if (!"untsq" %in% class(object)) {
    stop("'object' doesn't inherit class 'untsq'")
  } 
  invisible(object)
}

#'
validate_simsq <- function(object) {
  validate_sq(object)
  if (!"simsq" %in% class(object)) {
    stop("'object' doesn't inherit class 'simsq'")
  } 
  alph <- .get_alph(object)
  if (!all(alph %in% c(letters, "-"))) {
    stop("attribute 'alphabet' doesn't follow groups naming convention (lower latin letters and symbol '-')")
  }
  invisible(object)
}

#'
validate_atpsq <- function(object) {
  validate_sq(object)
  if (!"atpsq" %in% class(object)) {
    stop("'object' doesn't inherit class 'atpsq'")
  } 
  invisible(object)
}