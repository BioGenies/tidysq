#' @import tibble 
#' @exportClass sqtbl
#' @export
construct_sqtibble <- function(sq, name = NULL, type = NULL) {
  if (!is.null(type) &&
      !type %in% c("unt", "ami", "nuc")) {
    stop("'type' has to be one of 'unt', 'ami', 'nuc' or NULL")
  }
  if (is.list(sq)) {
    validate_sq(sq)
    sqtype <- .get_sq_type(sq)
    if (!is.null(type) && 
        type != sqtype) {
      stop("'type' given isn't identical to type 'sq' has; you don't need to pass 'type' argument")
    }
    validate_sq(sq, sqtype)
  } else if (is.character(sq)) {
    if (is.null(type)) {
      type <- "unt"
    }
    sq <- construct_sq(sq, type)
  } else {
    stop("'sq' should be either a vector of sequences or object of class 'sq'")
  }
  if (is.null(name)) {
    name <- names(sq)
  } else if (!is.character(name) ||
             !(length(name) == length(sq))){
    stop("'name' has to be a character vector of length equal to 'sq' or be NULL")
  }
  
  if (is.null(name)) {
    object <- tibble(sq = sq)
  } else {
    object <- tibble(name = name, sq = sq)
  }
  
  class(object) <- c("sqtbl", class(object))
  object
}

#'@import tibble 
validate_sqtibble <- function(object) {
  if (!"sqtbl" %in% class(object)) {
    stop("'object' doesn't inherit class 'sqtbl'")
  } 
  if (!"tbl" %in% class(object)) {
    stop("'object' doesn't inherit class 'tbl'")
  } 
  validate_tibble(object)
  sqcols_inds <- sapply(object, function(column) 'sq' %in% class(column))
  if (!any(sqcols_inds)) {
    stop("'object' contains no column of class 'sq'")
  }
  sqcols <- object[sqcols_inds]
  lapply(sqcols, validate_sq)
  invisible(object)
}
