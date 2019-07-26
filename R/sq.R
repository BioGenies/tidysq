#' @exportClass sq
#' @export
construct_sq <- function(sq, type = "unt", is_clean = NULL) {
  .check_sqstr_proper_char(sq)
  if (getOption("tidysq_no_check_mode") == TRUE) {
    .nc_construct_sq(sq, type, is_clean)
  } else {
    .check_type_in_ami_nuc_unt(type)
    .check_is_clean_in_TRUE_FALSE_NULL(is_clean)
    
    switch (type,
            ami = construct_amisq(sq),
            nuc = construct_nucsq(sq),
            unt = construct_untsq(sq))
  }
}

#' @export
.nc_construct_sq <- function(sq, type, is_clean) {
  .check_nc_type_in_ami_nuc(type)
  .check_nc_is_clean_in_TRUE_FALSE(is_clean)
  
  if (type == "ami") {
    if (is_clean) {
      sq <- .nc_bitify_sq_cami(sq)
      class(sq) <- c("clnsq", "amisq", "sq")
      attr(sq, "alphabet") <- aminoacids_df[!aminoacids_df[["amb"]], "one"]
      sq
    } else {
      sq <- .nc_bitify_sq_ami(sq)
      class(sq) <- c("amisq", "sq")
      attr(sq, "alphabet") <- aminoacids_df[["one"]]
      sq
    }
  } else if (type == "nuc") {
    if (is_clean) {
      sq <- .nc_bitify_sq_cnuc(sq)
      class(sq) <- c("clnsq", "nucsq", "sq")
      attr(sq, "alphabet") <- nucleotides_df[!nucleotides_df[["amb"]], "one"]
      sq
    } else {
      sq <- .nc_bitify_sq_nuc(sq)
      class(sq) <- c("nucsq", "sq")
      attr(sq, "alphabet") <- nucleotides_df[["one"]]
      sq
    }
  } 
}


#' @exportClass amisq
construct_amisq <- function(sq) {
  stop("not implemented!")
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  alph <- aminoacids_df[,"one"]
  is_ami_sq <- all(real_alph %in% alph)
  if (!is_ami_sq) {
    stop("each of letters should be in aminoacids alphabet (one of aminoacids_df[,'one'])")
  }
  
  object <- .bitify_sq(sq, alph)
  attr(object, "alphabet") <- alph
  class(object) <- c("amisq", "sq")
  object
}

#' @exportClass nucsq
construct_nucsq <- function(sq) {
  stop("not implemented!")
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  alph <- nucleotides_df[,"one"]
  is_nuc_sq <- all(real_alph %in% alph)
  if (!is_nuc_sq) {
    stop("each of letters should be in nucleotide alphabet (one of nucleotides_df[,'one'])")
  }
  
  object <- .bitify_sq(sq, alph)
  attr(object, "alphabet") <- alph
  class(object) <- c("nucsq", "sq")
  object
}

#' @exportClass untsq
construct_untsq <- function(sq) {
  stop("not implemented!")
  alph <- .get_real_alph(sq)
  
  object <- .bitify_sq(sq, alph)
  attr(object, "alphabet") <- alph
  class(object) <- c("untsq", "sq")
  object
}

validate_sq <- function(object, type = NULL) {
  if (!"sq" %in% class(object))
    stop("'object' doesn't inherit class 'sq'")
  sqtype <- .get_sq_subclass(object)
  if (!(length(sqtype) == 1)) 
    stop("'object' should have exactly one of types: 'ami', 'nuc', 'unt', 'sim', 'atp")
  if (is.null(attr(object, "alphabet"))) 
    stop("'object' doesn't 'alphabet' attribute")
  if (is.null(attr(object, "alphabet")))
    stop("'object' doesn't 'alphabet' attribute")
  alph <- .get_alph(object)
  if (!is.character(alph) &&
      !is.numeric(alph))
    stop("attribute 'alphabet' is neither a character nor a numeric vector")
  #assumption about length of one of each character - this can be changed in future
  if (!all(sapply(alph, length) == 1))
    stop("attribute 'alphabet' have elements that aren't one element long")
  if (!is.list(object))
    stop("'object' isn't a list")
  if (!all(sapply(object, is.raw))) 
    stop("'object' isn't a list of raw vectors")
  # if (???) {
  #   stop("'alphabet' attribute has less elements than are different values in 'object'")
  # } - quite long step, is it necessary?
  if (!is.null(type)) {
    switch (type,
      ami = validate_amisq(object),
      nuc = validate_nucsq(object),
      unt = validate_untsq(object),
      sim = validate_simsq(object),
      atp = validate_atpsq(object),
      enc = validate_encsq(object)
    )
  }
  
  invisible(object)
}

validate_nucsq <- function(object) {
  if (!"nucsq" %in% class(object))
    stop("'object' doesn't inherit class 'nucsq'")
  alph <- .get_alph(object)
  if (!("clnsq" %in% class(object))) {
    if (!identical(alph, nucleotides_df[,"one"]))
      stop("attribute 'alphabet' isn't identical to standard nucleotides alphabet")
  } else if (!identical(alph, nucleotides_df[!nucleotides_df[["amb"]],"one"]))
      stop("attribute 'alphabet' isn't identical to cleaned nucleotides alphabet")
  
  invisible(object)
}

validate_amisq <- function(object) {
  if (!"amisq" %in% class(object)) 
    stop("'object' doesn't inherit class 'amisq'")
  alph <- .get_alph(object)
  if (!("clnsq" %in% class(object))) {
    if (!identical(alph, aminoacids_df[,"one"]))
      stop("attribute 'alphabet' isn't identical to standard aminoacids alphabet")
  } else if (!identical(alph, aminoacids_df[!aminoacids_df[["amb"]],"one"])) 
      stop("attribute 'alphabet' isn't identical to cleaned aminoacids alphabet")
  
  
  invisible(object)
}

validate_untsq <- function(object) {
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector")
  if (!"untsq" %in% class(object))
    stop("'object' doesn't inherit class 'untsq'")
  
  invisible(object)
}

validate_simsq <- function(object) {
  if (!"simsq" %in% class(object))
    stop("'object' doesn't inherit class 'simsq'")
  alph <- .get_alph(object)
  if (!all(alph %in% c(letters, "-"))) 
    stop("attribute 'alphabet' doesn't follow groups naming convention (lower latin letters and symbol '-')")
  
  invisible(object)
}

validate_atpsq <- function(object) {
  if (!"atpsq" %in% class(object))
    stop("'object' doesn't inherit class 'atpsq'")
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector")
  
  invisible(object)
}

validate_encsq <- function(object) {
  if (!"encsq" %in% class(object))
    stop("'object' doesn't inherit class 'encsq'")
  alph <- .get_alph(object)
  if (!is.numeric(alph))
    stop("attribute 'alphabet' isn't a numeric vector")
  
  invisible(object)
}