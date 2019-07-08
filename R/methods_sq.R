#methods other than print of class sq

#' @exportMethod `[` sq
#' @export
`[.sq` <- function(i, ...) {
  ret <- NextMethod()
  class(ret) <- "sq"
  ret
}