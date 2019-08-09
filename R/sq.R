#' @exportClass sq
#' @export
construct_sq <- function(sq, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_sqstr_proper_char(sq)
  if (getOption("tidysq_no_check_mode") == TRUE) {
    .nc_construct_sq(sq, type, is_clean)
  } else {
    .check_type_or_nonst_alph(type, is_clean, non_standard)
    if (!is.null(non_standard)) {
      .nonst_construct_sq(sq, non_standard)
    } else {
      .check_is_clean_in_TRUE_FALSE_NULL(is_clean)
      
      if (is.null(type)) {
        type <- .guess_sq_type(sq)
      }
      .check_type_in_ami_nuc_unt(type)
      
      switch (type,
              ami = construct_amisq(sq, is_clean),
              nuc = construct_nucsq(sq, is_clean),
              unt = construct_untsq(sq))
    }
  }
}

.nc_construct_sq <- function(sq, type, is_clean) {
  .check_nc_type_in_ami_nuc(type)
  .check_nc_is_clean_in_TRUE_FALSE(is_clean)
  
  sq <- .nc_bitify_sq(sq, type, is_clean)
  sq <- .set_class(sq, type, is_clean)
  .set_alph(sq, .get_standard_alph(type, is_clean))
}


#' @importFrom stringi stri_sub
#' @importFrom stringi stri_locate_all_regex
.nonst_construct_sq <- function(sq, non_standard) {
  .check_nonst_proper_char(non_standard)
  .check_nonst_nchar(non_standard)
  
  sq <- lapply(sq, function(s) {
    pos <- stri_locate_all_regex(s, non_standard)
    binded <- apply(na.omit(do.call(rbind, pos)), 2, sort)
    if (length(binded) == 0) {
      strsplit(s, "")[[1]]
    } else {
      if (length(binded) == 2) {
        binded <- t(as.matrix(binded))
        res_ind <- binded[1]:binded[2]
      } else res_ind <- as.integer(binded, 1, function(row) row[1]:row[2])
      sin_ind <- setdiff(1:nchar(s), res_ind)
      ind <- .merge_ind(sin_ind, binded[,1])
      n <- nrow(binded) + length(sin_ind)
      begs <- integer(n)
      ends <- integer(n)
      begs[(1:n)[ind]] <- sin_ind
      ends[(1:n)[ind]] <- sin_ind
      begs[(1:n)[!ind]] <- binded[,1]
      ends[(1:n)[!ind]] <- binded[,2]
      stri_sub(s, begs, ends)
    }
  })
  
  alph <- unique(unlist(sq))
  .check_alph_length(alph)
  sq <- .bitify_sq(sq, alph)
  sq <- .set_alph(sq, alph)
  .set_class(sq, "atp")
}

#' @exportClass amisq
construct_amisq <- function(sq, is_clean) {
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  if (!is.null(is_clean) &&
      is_clean == TRUE &&
      !all(real_alph %in% .get_standard_alph("ami", TRUE))) {
    stop("'is_clean' is given TRUE, but sequences contain at least one ambiguous aminoacid")
  }
  if (is.null(is_clean)) {
    is_clean <- .guess_ami_is_clean(real_alph)
  }
  sq <- .nc_bitify_sq(sq, "ami", is_clean)
  sq <- .set_alph(sq, .get_standard_alph("ami", is_clean))
  .set_class(sq, "ami", is_clean)
}

#' @exportClass nucsq
construct_nucsq <- function(sq, is_clean) {
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  if (!is.null(is_clean) &&
      is_clean == TRUE &&
      !all(real_alph %in% .get_standard_alph("nuc", TRUE))) {
    stop("'is_clean' is given TRUE, but sequences contain at least one ambiguous nucleotide")
  }
  if (is.null(is_clean)) {
    is_clean <- .guess_nuc_is_clean(real_alph)
  }
  sq <- .nc_bitify_sq(sq, "nuc", is_clean)
  sq <- .set_alph(sq, .get_standard_alph("nuc", is_clean))
  .set_class(sq, "nuc", is_clean)
}

#' @exportClass untsq
construct_untsq <- function(sq) {
  alph <- .get_real_alph(sq)
  .check_alph_length(alph)
  
  sq <- .bitify_sq(sq, alph)
  sq <- .set_alph(sq, alph)
  .set_class(sq, "unt", FALSE)
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
    if (!identical(alph, .get_standard_alph("nuc", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard nucleotides alphabet")
  } else if (!identical(alph, .get_standard_alph("nuc", TRUE)))
      stop("attribute 'alphabet' isn't identical to cleaned nucleotides alphabet")
  
  invisible(object)
}

validate_amisq <- function(object) {
  if (!"amisq" %in% class(object)) 
    stop("'object' doesn't inherit class 'amisq'")
  alph <- .get_alph(object)
  if (!("clnsq" %in% class(object))) {
    if (!identical(alph, .get_standard_alph("ami", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard aminoacids alphabet")
  } else if (!identical(alph, .get_standard_alph("ami", TRUE))) 
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