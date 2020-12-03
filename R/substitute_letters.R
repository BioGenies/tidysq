#' Substitute letters in a sequence
#' 
#' @description 1) Replace ambiguous/extraordinary letters in a nucleic or 
#' amino acid sequence, stored in a \code{\link[=sq-class]{sq}} object, with the ones
#' that are compliant with the IUPAC standard, ones that are user-defined 
#' or with \code{NA} values.
#' 
#' 2) Replace default amino acid letters in a sequence with a custom encoding 
#' to create sequences with simplified alphabets.
#' 
#' The function is only used to replace letters in the alphabet. 
#' It cannot be used to merge multiple characters into one.
#' 
#' @template x
#' @param encoding [\code{character}]\cr
#'  Letters to be replaced together with their replacements.
#'  One letter can be replaced with multiple symbols.
#'  To perform substitution create a named vector, e.g.
#'  \code{c(A = Ala, H = His, amino_or_nucleic_acid_symbol = replacement)}.
#' @template NA_letter
#' @template three-dots
#' 
#' @return a \code{\link[=sq-class]{sq}} object with \strong{atp} type with replaced alphabet,
#' defined by user.
#' 
#' @details \code{substitute_letters} allows to replace ambiguous/extraordinary 
#' letters in nucleic or amino acid sequence with user-defined or IUPAC 
#' symbols. Letters can also be replaced with \code{\link{NA}} values, so that they 
#' can be later removed, from the sequence, by \code{\link{clean}} function.
#' 
#' \code{substitute_letters} can be used to replace default amino acid letters 
#' with encodings. They can be user-defined or be derived from various 
#' simplified alphabets.
#' 
#' One letter of the alphabet may be replaced by a multiple character. 
#' 
#' The alphabet characters to be replaced need to be written in capital letters
#' and must originate from default alphabets, otherwise error will be 
#' introduced.
#' 
#' Multiple string of letters to be substituted 
#' (ex. \code{c(AHG = "replacement")}) will also produce an error.
#' 
#' Replacing multiple letters with the same symbol 
#' (ex. \code{c(A = "rep1", H  = "rep1", G = "rep1")}) is allowed.
#' 
#' Created sequence will be stripped of \strong{cln} subtype, 
#' if the original sequence possessed it. This will also occur when
#' the letter to be replaced will not be found in the sequence. 
#' It remain unchanged but will lose subclass.
#' 
#' The newly constructed will have a new type \strong{atp}, 
#' representing sequences with atypical alphabet.
#' 
#' All replaced letters will have the character type. 
#' Multiple symbol replacement will be recognized as one letter and 
#' cannot be separated in future operations into single letters. 
#' 
#' @seealso \code{\link[=sq-class]{sq}}
#' 
#' @export
substitute_letters <- function(x, encoding, ...)
  UseMethod("substitute_letters")

#' @export
substitute_letters.default <- function(x, encoding, ...)
  stop("cannot substitute letters in this type of object", call. = FALSE)

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
