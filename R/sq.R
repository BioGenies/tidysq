#' @exportClass nucsq
#' @export
construct_nucsq <- function(sq) {
  if (!is.character(sq)) {
    stop("'sq' has to be a character vector", call. = FALSE)
  }
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    if (!all(nchar(sq) == 1)) {
      stop("'sq' should have one element with whole sequence or many elements of length 1", call. = FALSE)
    }
  }
  is_nuc_sq <- all(sq %in% nucleotides_df[,"one"])
  if (!is_nuc_sq) {
    stop("each of sequence elements should be in nucleotide alphabet (one of nucleotides_df['one'])", call. = FALSE)
  }
  
  object <- factor(sq, levels = nucleotides_df[,"one"])
  class(object) <- c("nucsq", "sq", class(object))
  object
}

#' @exportClass aasq
#' @export
construct_aasq <- function(sq) {
  # TO DO: what if user gives list of aminoacids three letter names
  if (!is.character(sq)) {
    stop("'sq' has to be a character vector", call. = FALSE)
  }
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    if (!all(nchar(sq) == 1)) {
      stop("'sq' should have one element with whole sequence or many elements of length 1", call. = FALSE)
    }
  }
  is_aa_sq <- all(sq %in% aminoacids_df[,"one"])
  if (!is_aa_sq) {
    stop("each of sequence elements should be in aminoacids alphabet (latin letters with - or .)", call. = FALSE)
  }
  
  object <- factor(sq, levels = aminoacids_df[,"one"])
  class(object) <- c("aasq", "sq", class(object))
  object
}

#' @exportClass untsq
#' @export
construct_untsq <- function(sq) {
  if (!is.character(sq)) {
    stop("'sq' has to be a character vector", call. = FALSE)
  }
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    if (!all(nchar(sq) == 1)) {
      stop("'sq' should have one element with whole sequence or many elements of length 1", call. = FALSE)
    }
  }
  
  object <- factor(sq, levels = sort(unique(sq)))
  class(object) <- c("untsq", "sq", class(object))
  object
}

#' @exportClass nucsq
#' @export
construct_sq <- function(sq, type = "unt") {
  if (!type %in% c("unt", "aa", "nuc")) {
    stop("'type' has to be one of 'unt', 'aa', 'nuc'", call. = FALSE)
  }
  
  if (type == "aa") {
    construct_aasq(sq)
  } else if (type =="nuc") {
    construct_nucsq(sq)
  } else {
    construct_untsq(sq)
  }
}

#'
validate_sq <- function(object) {
  if (!"sq" %in% class(object)) {
    stop("'object' doesn't inherit class 'sq'")
  } 
  if (!"factor" %in% class(object)) {
    stop("'object' doesn't inherit class 'factor'")
  }
  invisible(object)
}

#'
validate_aasq <- function(object) {
  validate_sq(object)
  if (!"aasq" %in% class(object)) {
    stop("'object' doesn't inherit class 'aasq'")
  } 
  if (!all(levels(object) == aminoacids_df[,"one"])) {
    stop("'object' levels aren't identical to standard aminoacids alphabet")
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
validate_nucsq <- function(object) {
  validate_sq(object)
  if (!"nucsq" %in% class(object)) {
    stop("'object' doesn't inherit class 'nucsq'")
  } 
  if (!all(levels(object) == nucleotides_df[,"one"])) {
    stop("'object' levels aren't identical to standard nucleotides alphabet")
  }
  invisible(object)
}