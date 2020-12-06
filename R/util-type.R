#' Get type of an sq object
#' 
#' @description Returns type of sequences/alphabet contained in
#' \code{\link[=sq-class]{sq}} object.
#'  
#' @template x
#' @template three-dots
#'  
#' @return A string, one of: "ami_bsc", "ami_ext", "dna_bsc", "dna_ext",
#' "rna_bsc", "rna_ext", "unt" or "atp".
#' 
#' @details
#' Types returned by this function can be passed as argument to functions like
#' \code{\link{random_sq}} and \code{\link{find_invalid_letters}}.
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGA", "TGACGAGCTTA", "ACTTTAGC"),
#'              alphabet = "dna_bsc")
#'
#' # Extracting type of sq objects:
#' sq_type(sq_ami)
#' sq_type(sq_dna)
#'
#' # Classes are tightly related to these types:
#' class(sq_ami)[1]
#' class(sq_dna)[1]
#'
#' @family type_functions
#' @seealso \code{\link[=sq-class]{sq class}}
#' @export
sq_type <- function(x, ...)
  UseMethod("sq_type")

#' @export
sq_type.default <- function(x, ...)
  stop("cannot determine sq_type of this type of object", call. = FALSE)

#' @rdname sq_type
#' @export
sq_type.sq <- function(x, ...)
  vec_ptype_abbr(x)

#' @rdname sq_type
#' @export
`sq_type<-` <- function(x, value)
  UseMethod("sq_type<-")

#' @export
`sq_type<-.default` <- function(x, value)
  stop("cannot change sq_type of this type of object", call. = FALSE)

#' @rdname sq_type
#' @export
`sq_type<-.sq` <- function(x, value) {
  typify(x, value)
}

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
