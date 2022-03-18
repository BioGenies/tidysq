#' Substitute letters in a sequence
#' 
#' @description Replaces all occurrences of a letter with another.
#' 
#' @template x
#' @param encoding [\code{character} || \code{numeric}]\cr
#'  A dictionary (named vector), where names are letters to be replaced and
#'  elements are their respective replacements.
#' @template NA_letter
#' @template three-dots
#' 
#' @return An \code{\link[=sq-class]{sq}} object of \strong{atp} type with
#' updated alphabet.
#' 
#' @details
#' \code{substitute_letters} allows to replace unwanted letters in any sequence
#' with user-defined or IUPAC  symbols. Letters can also be replaced with
#' \code{\link{NA}} values, so that they  can be later removed from the sequence
#' by \code{\link{remove_na}} function.
#'
#' It doesn't matter whether replaced or replacing letter is single or multiple
#' character. However, the user cannot replace multiple letters with one nor one
#' letter with more than one.
#'
#' Of course, multiple different letters can be encoded to the same symbol, so
#' \code{c(A = "rep1", H = "rep1", G = "rep1")} is allowed, but
#' \code{c(AHG = "rep1")} is not (unless there is a letter "\code{AHG}" in
#' the alphabet). By doing that any information of separateness of original
#' letters is lost, so it isn't possible to retrieve original sequence after
#' this operation.
#'
#' All encoding names must be letters contained within the alphabet, otherwise
#' an error will be thrown.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
#'              alphabet = "dna_ext")
#' sq_ami <- sq(c("MIOONYTWIL","TIOOLGNIIYROIE", "NYERTGHLI", "MOYXXXIOLN"),
#'              alphabet = "ami_ext")
#' sq_atp <- sq(c("mALPVQAmAmA", "mAmAPQ"), alphabet = c("mA", LETTERS))
#'
#' # Not all letters must have their encoding specified:
#' substitute_letters(sq_dna, c(T = "t", A = "a", C = "c", G = "g"))
#' substitute_letters(sq_ami, c(M = "X"))
#'
#' # Multiple character letters are supported in encodings:
#' substitute_letters(sq_atp, c(mA = "-"))
#' substitute_letters(sq_ami, c(I = "ough", O = "eau"))
#'
#' # Numeric substitutions are allowed too, these are coerced to characters:
#' substitute_letters(sq_dna, c(N = 9, G = 7))
#'
#' # It's possible to replace a letter with NA value:
#' substitute_letters(sq_ami, c(X = NA_character_))
#'
#' @family type_functions
#' @export
substitute_letters <- function(x, encoding, ...)
  UseMethod("substitute_letters")

#' @export
substitute_letters.default <- function(x, encoding, ...)
  stop_no_method(substitute_letters, x)

#' @rdname substitute_letters
#' @export
substitute_letters.sq <- function(x, encoding, ...,
                                  NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  assert_atomic_vector(encoding, names = "unique")
  assert_subset(names(encoding), alphabet(x))
  
  if (is.numeric(encoding)) {
    assert_integerish(encoding)
    # Changes storage mode, because it preserves attributes
    # (and we want to preserve names)
    mode(encoding) <- "character"
  } else if (is.character(encoding)) {
    assert_character(encoding)
  } else {
    stop("encoding must be either numeric of character vector", call. = FALSE)
  }
  
  CPP_substitute_letters(x, encoding, NA_letter)
}
