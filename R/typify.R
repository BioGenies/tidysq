#' Set type of an sq object
#'
#' @description Sets sequence type (and, consequently, alphabet attribute) to
#' one of \strong{ami}, \strong{dna} or \strong{rna} types.
#'
#' @template x
#' @template dest_type
#' @template NA_letter
#' @template three-dots
#' 
#' @return \code{\link[=sq-class]{sq}} object with the same letters as input
#' \code{x}, but with type as specified in \code{dest_type}.
#' 
#' @details
#' Sometimes functions from I/O module return sequences of incorrect type, most
#' often \strong{unt} (which indicates no type). It happens mostly whenever
#' there are letters that don't fit into target alphabet. After replacing wrong
#' letters with correct ones with \code{\link{substitute_letters}} the user has
#' sequences of type \strong{atp}, even if their alphabet is contained in the
#' target one. At the same time, many functions demand sequences to be of
#' standard type (i.e. \strong{ami}, \strong{dna} or \strong{rna}) or behave
#' differently for these.
#'
#' \code{typify()} is used to help with these situations by allowing the user
#' to convert their sequences to target type. There are some conditions that
#' must be met to use this function. The most important is that typified
#' \code{sq} object must not contain invalid letters. If this condition is not
#' satisfied, an error is thrown.
#'
#' If \code{dest_type} is equal to type of \code{sq}, function simply returns
#' input value.
#'
#' @examples
#' # Constructing sq object with strange characters (type will be set to "unt"):
#' sq_unt <- sq(c("&VPLG&#", "##LCG"))
#'
#' # Substituting letters with "X", which stands for unknown amino acid:
#' sq_sub <- substitute_letters(sq_unt, c(`&` = "X", `#` = "X"))
#'
#' # Setting extended amino acid type (only extended one has "X" letter):
#' typify(sq_sub, "ami_ext")
#'
#' @family type_functions
#' @export
typify <- function(x, dest_type, ...) {
  assert_sq_type(dest_type)
  
  UseMethod("typify")
}

#' @export
typify.default <- function(x, dest_type, ...)
  stop("'typify' isn't implemented for this type of object", call. = FALSE)

#' @rdname typify
#' @export
typify.sq <- function(x, dest_type, ...,
                      NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_typify(x, dest_type, NA_letter)
}
