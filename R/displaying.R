#' @importFrom pillar type_sum
#' @exportMethod type_sum sqcol
#' @export
type_sum.sqcol <- function(x) {
  "sq"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum sq
#' @export
type_sum.sq <- function(x) {
  "sq"
}

#' @importFrom pillar is_vector_s3
#' @exportMethod is_vector_s3 sqcol
#' @export
is_vector_s3.sqcol <- function(x) {
  TRUE
}

obj_sum.sq <- function(x) {
  "squauauau"
}

#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft_simple
#' @exportMethod pillar_shaft sqcol
#' @export
pillar_shaft.sqcol <- function(x, ...) {
  out <- paste0("<", as.character(sapply(x, function(i) as.character(i[[1]]))), ">")
  new_pillar_shaft_simple(out, align = "right")
}
#'