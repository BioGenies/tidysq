#' Get type of an sq object
#' 
#' @description Returns type of sequences/alphabet contained in
#' \code{\link[=sq-class]{sq}} object.
#'  
#' @template x
#' @template three-dots
#' @param value [\code{character(1)}]\cr
#'  The name of destination type - any valid sq type.
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
  stop_no_method(
    sq_type, x,
    msg = function(cls) paste0("cannot determine sq type of object of classes <",
                               paste0(class(x), collapse = ", "), ">")
  )

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
  stop_no_method(
    `sq_type<-`, x,
    msg = function(cls) paste0("cannot change sq type of object of classes <",
                               paste0(class(x), collapse = ", "), ">")
  )

#' @rdname sq_type
#' @export
`sq_type<-.sq` <- function(x, value) {
  typify(x, value)
}

interpret_type <- function(name) {
  name <- tolower(gsub(" ", "_", name))
  if      (name %in% c("dna_bsc", "bsc_dna", "basic_dna", "dna_basic")) "dna_bsc"
  else if (name %in% c("dna_ext", "ext_dna", "extended_dna", "dna_extended", "dna")) "dna_ext"
  else if (name %in% c("rna_bsc", "bsc_rna", "basic_rna", "rna_basic")) "rna_bsc"
  else if (name %in% c("rna_ext", "ext_rna", "extended_rna", "rna_extended", "rna")) "rna_ext"
  else if (name %in% c("ami_bsc", "bsc_ami", "basic_ami", "ami_basic")) "ami_bsc"
  else if (name %in% c("ami_ext", "ext_ami", "extended_ami", "ami_extended", "ami", "aminoacids")) "ami_ext"
  else if (name %in% c("unt", "untyped")) "unt"
  else if (name %in% c("atp", "atypical")) stop("When creating atp sq, alphabet should be vector of letters", call. = FALSE)
  else stop("Cannot interpret type for provided alphabet", call. = FALSE)
}

type_as_class <- function(type)
  paste0("sq_", type)
