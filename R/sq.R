#' @exportClass sq
#' @export
construct_sq <- function(sq, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_character(sq, "'sq'", allow_zero_len = TRUE)
  if (getOption("tidysq_no_check_mode") == TRUE) {
    .nc_construct_sq(sq, type, is_clean)
  } else {
    .check_type_or_nonst_alph(type, is_clean, non_standard)
    if (!is.null(non_standard)) {
      .nonst_construct_sq(sq, non_standard)
    } else {
      .check_logical(is_clean, "'is_clean'", allow_null = TRUE, single_elem = TRUE)
      if (is.null(type)) {
        type <- .guess_sq_type(sq)
      }
      .check_type(type, allow_unt = TRUE)
      switch(type,
             ami = .construct_amisq(sq, is_clean),
             nuc = .construct_nucsq(sq, is_clean),
             unt = .construct_untsq(sq))
    }
  }
}

.nc_construct_sq <- function(sq, type, is_clean) {
  .check_type(type)
  .check_logical(is_clean, single_elem = TRUE)
  
  sq <- .nc_bitify_sq(sq, type, is_clean)
  sq <- .set_class(sq, type, is_clean)
  .set_alph(sq, .get_standard_alph(type, is_clean))
}


#' @importFrom stringi stri_sub
#' @importFrom stringi stri_locate_all_regex
.nonst_construct_sq <- function(sq, non_standard) {
  .check_character(non_standard, "'non_standard'")
  .check_nchar(non_standard, "'non_standard'", minimal_nchar = 2)
  
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
.construct_amisq <- function(sq, is_clean) {
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  .check_real_alph_clean(real_alph, "ami", is_clean)
  if (is.null(is_clean)) {
    is_clean <- .guess_ami_is_clean(real_alph)
  }
  sq <- .nc_bitify_sq(sq, "ami", is_clean)
  sq <- .set_alph(sq, .get_standard_alph("ami", is_clean))
  .set_class(sq, "ami", is_clean)
}

#' @exportClass nucsq
.construct_nucsq <- function(sq, is_clean) {
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  .check_real_alph_clean(real_alph, "nuc", is_clean)
  if (is.null(is_clean)) {
    is_clean <- .guess_nuc_is_clean(real_alph)
  }
  sq <- .nc_bitify_sq(sq, "nuc", is_clean)
  sq <- .set_alph(sq, .get_standard_alph("nuc", is_clean))
  .set_class(sq, "nuc", is_clean)
}

#' @exportClass untsq
.construct_untsq <- function(sq) {
  alph <- .get_real_alph(sq)
  .check_alph_length(alph)
  
  sq <- .bitify_sq(sq, alph)
  sq <- .set_alph(sq, alph)
  .set_class(sq, "unt", FALSE)
}

validate_sq <- function(object, type = NULL) {
  argname <- deparse(substitute(object))
  if (!"sq" %in% class(object))
    stop(argname, " doesn't inherit class 'sq'", call. = FALSE)
  sqtype <- .get_sq_subclass(object)
  if (length(sqtype) != 1) 
    stop("'object' should have exactly one of classes: 'amisq', 'nucsq', 'untsq', 'encsq', 'atpsq'", call. = FALSE)
  if (is.null(attr(object, "alphabet"))) 
    stop("'object' doesn't 'alphabet' attribute", call. = FALSE)
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
           unt = .validate_untsq(object),
           atp = .validate_atpsq(object),
           enc = .validate_encsq(object),
           invisible(object)
    )
  }
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
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector", call. = FALSE)
  if (!"untsq" %in% class(object))
    stop("'object' doesn't inherit class 'untsq'", call. = FALSE)
  
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