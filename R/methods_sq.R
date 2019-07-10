#methods other than print of class sq

#' @exportMethod `[` sq
#' @export
`[.sq` <- function(x, i, j, ...) {
  ret <- NextMethod()
  class(ret) <- class(x)
  attr(ret, "alphabet") <- .get_alph(x)
  ret
}

#' @exportMethod as.character sq
#' @export
as.character.sq <- function(x, ...) {
  dict <- .get_alph(x)
  sapply(x, function(s) paste(ifelse(!is.na(dict[s]), dict[s], "*"), collapse = ""))
}