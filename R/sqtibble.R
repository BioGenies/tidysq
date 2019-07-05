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
  sqtypes <- sapply(sq, function(s) intersect(class(s), c("aasq", "nucsq", "simsq", "untsq")))
  sqtypes <- unlist(sqtypes)
  if (length(setdiff(sqtypes, "untsq")) > 1) {
    handle_opt_txt("tidysq_constr_mtype_action",
                   "column 'sq' contains more than one type of sequences among 'nuc', 'aa' and 'sim' - not all functions can handle that")
  }
  if (any(sqtypes == "untsq")) {
    handle_opt_txt("tidysq_constr_unt_action",
                   "column 'sq' contains at least one 'unt' sequence - not all functions can handle that")
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
  if (!is.character(object[["name"]])) {
    stop("column 'name' is not a character column")
  }
  if (!is.list(object[["sq"]])) {
    stop("column 'sq' is not a list")
  }
  if (class(object[["sq"]]) != "sqcol") {
    stop("column 'sq' is not kept as 'sqcol'")
  }
  sapply(object[["sq"]], validate_sq)
  invisible(object)
}
