#' @import tibble::tiblle
#' @exportClass sqtibble
#' @export
construct_sqtibble <- function(name, sq) {
  if (!is.character(name)) {
    stop("name should be character vector")
  }
  # after adding validate :
  # if (!(is.list(sq) &&
  #       all(validate_sq(sq)))) {
  #   stop("sq should be list of sq objects")
  # }
  if (any(isClass("aasq")) && any(isClass("nucsq"))) {
    #later it could be a option of package:
    warning("sq contains both nuc and aa sequences - not all functions support that")
  }
  if (any(isClass("ambsq"))) {
    #later it could be a option of package:
    warning("sq contains at least one amb sequence - not all functions support that")
  }
  object <- tibble(name = name, sq = sq)
  class(object) <- c("sqtibble", class(object))
  object
}
