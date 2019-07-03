#' @exportClass nucsq
#' @export
construct_nucsq <- function(sq) {
  if (!is.character(sq)) {
    stop("sq must be character vector", call. = FALSE)
  }
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    if (all(nchar(sq) == 1)) {
      stop("sq should have one element with whole sequence or many elements of length 1", call. = FALSE)
    }
  }
  is_nuc_sq <- all(sq %in% nucleotides_df[,"one"])
  if (!is_nuc_sq) {
    stop("each of sequence elements should be in standard nucleotide alphabet", call. = FALSE)
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
    stop("sq must be character vector", call. = FALSE)
  }
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    if (all(nchar(sq) == 1)) {
      stop("sq should have one element with whole sequence or many elements of length 1", call. = FALSE)
    }
  }
  is_aa_sq <- all(sq %in% aminoacids_df[,"one"])
  if (!is_aa_sq) {
    stop("each of sequence elements should be in standard aminoacids alphabet", call. = FALSE)
  }
  
  object <- factor(sq, levels = aminoacids_df[,"one"])
  class(object) <- c("aasq", "sq", class(object))
  object
}

#' @exportClass ambsq
#' @export
construct_ambsq <- function(sq) {
  if (!is.character(sq)) {
    stop("sq must be character vector", call. = FALSE)
  }
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    if (all(nchar(sq) == 1)) {
      stop("sq should have one element with whole sequence or many elements of length 1", call. = FALSE)
    }
  }
  #should this be all?
  is_amb_sq <- all(sq %in% c(LETTERS, "-", "."))
  if (!is_amb_sq) {
    stop("each of sequence elements should be in extended alphabet (latin letters or +, .)", call. = FALSE)
  }
  
  object <- factor(sq, levels = c(LETTERS, "-", "."))
  class(object) <- c("ambsq", "sq", class(object))
  object
}

#' @exportClass nucsq
#' @export
construct_sq <- function(sq, type = "amb") {
  if (!type %in% c("amb", "aa", "nuc")) {
    stop("type must be one of 'amb', 'aa', 'nuc'", call. = FALSE)
  }
  
  if (type == "aa") {
    construct_aasq(sq)
  } else if (type =="nuc") {
    construct_nucsq(sq)
  } else {
    construct_ambsq(sq)
  }
}

#'
validate_sq <- function(object) {
  if (!"sq" %in% class(object)) {
    stop("object doesn't inherit class 'sq'")
  } 
  if (!"factor" %in% class(object)) {
    stop("object doesn't inherit class 'factor'")
  }
  if (!all(levels(object) %in% c(LETTERS, "-", "."))) {
    stop("object elements should be elements of extended alphabet (latin letters or +, .)")
  }
  invisible(object)
}

#'
validate_aasq <- function(object) {
  validate_sq(object)
  if (!"aasq" %in% class(object)) {
    stop("object doesn't inherit class 'aasq'")
  } 
  if (!all(levels(object) == aminoacids_df[,"one"])) {
    stop("object elements should be identical to standard aminoacids alphabet")
  }
  invisible(object)
}

#'
validate_ambsq <- function(object) {
  validate_sq(object)
  if (!"ambsq" %in% class(object)) {
    stop("object doesn't inherit class 'ambsq'")
  } 
  if (!all(levels(object) %in% c(LETTERS, "-", "."))) {
    stop("object elements should be identical to extended alphabet (latin letters or +, .)")
  }
  invisible(object)
}

#'
validate_nucsq <- function(object) {
  validate_sq(object)
  if (!"nucsq" %in% class(object)) {
    stop("object doesn't inherit class 'nucsq'")
  } 
  if (!all(levels(object) == nucleotides_df[,"one"])) {
    stop("object elements should be identical to standard aminoacids alphabet")
  }
  invisible(object)
}