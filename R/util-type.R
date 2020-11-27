interpret_type <- function(name) {
  # TODO: improve; just improve
  switch(name, 
         dna_bsc = "dna_bsc",
         dna_ext = "dna_ext",
         rna_bsc = "rna_bsc",
         rna_ext = "rna_ext",
         ami_bsc = "ami_bsc",
         ami_ext = "ami_ext",
         unt = "unt")
}

type_as_class <- function(type)
  paste0("sq_", type)

#' Get type of a sq object
#' 
#' Function checks which type of sequences are contained in \code{\link[=sq-class]{sq}} object.
#'  
#' @param x a \code{\link[=sq-class]{sq}} object to be checked.
#' @template three-dots
#'  
#' @return A \code{\link{character}} string, type of\code{\link[=sq-class]{sq}} object - can be one of
#' "ami", "dna", "rna", "unt", "atp" or "enc".
#' 
#' @details This function returns type of sequence from \code{\link[=sq-class]{sq}} object.
#' If the type of sequence is \strong{dna}, \strong{rna}, \strong{ami}, \strong{unt},
#' \strong{atp} or \strong{enc} function returns "dna", "rna", "ami", "unt", "atp" or
#' "enc" respectivetly.
#'  
#' @seealso \code{\link[=sq-class]{sq}} \code{\link[=sq-class]{sq}}
#' @export
sq_type <- function(x, ...)
  UseMethod("sq_type")

#' @rdname sq_type
#' @export
sq_type.default <- function(x, ...)
  stop("cannot determine sq_type of this type of object", call. = FALSE)

#' @rdname sq_type
#' @export
sq_type.sq <- function(x, ...)
  vec_ptype_abbr(x)

#' @export
`sq_type<-` <- function(x, value)
  UseMethod("sq_type<-")

#' @export
`sq_type<-.default` <- function(x, value)
  stop("cannot change sq_type of this type of object", call. = FALSE)

#' @export
`sq_type<-.sq` <- function(x, value) {
  typify(x, value)
}
