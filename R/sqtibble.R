#' @import tibble 
#' @exportClass sqtbl
#' @exportClass sqcol
#' @export
construct_sqtibble <- function(name, sq) {
  if (!is.character(name)) {
    stop("name should be character vector")
  }
  if (!(is.list(sq))) {
    stop("sq should be list of 'sq' objects")
  }
  sapply(sq, validate_sq)
  class(sq) <- c("sqcol")
  if (any(sapply(sq, function(obj) "aasq" %in% class(obj))) && 
      any(sapply(sq, function(obj) "nucsq" %in% class(obj)))) {
    #later it could be a option of package:
    warning("column 'sq' contains both 'nuc' and 'aa' types sequences - not all functions can handle that", call. = FALSE)
  }
  if (any(sapply(sq, function(obj) "ambsq" %in% class(obj)))) {
    #later it could be a option of package:
    warning("column 'sq' contains at least one 'unt' sequence - not all functions can handle that", call. = FALSE)
  }
  object <- tibble(name = name, sq = sq)
  class(object) <- c("sqtbl", class(object))
  object
}

#'@import tibble 
validate_sqtibble <- function(object) {
  if (!"sqtbl" %in% class(object)) {
    stop("object doesn't inherit class 'sqtbl'")
  } 
  if (!"tbl" %in% class(object)) {
    stop("object doesn't inherit class 'tbl'")
  } 
  validate_tibble(object)
  if (!has_name(object, "name")) {
    stop("there's no column named 'name'")
  }
  if (!has_name(object, "sq")) {
    stop("there's no column named 'sq'")
  }
  if (!is.character(object$name)) {
    stop("column 'name' is not a character column")
  }
  if (!is.list(object$sq)) {
    stop("column 'sq' is not a list")
  }
  if (!is.list(object$sq)) {
    stop("column 'sq' is not kept as 'sqvec'")
  }
  sapply(object$sq, validate_sq)
  invisible(object)
}
