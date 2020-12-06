#' Remove sequences that contain NA values
#'
#' @description This function replaces sequences with \code{NA} values by
#' empty (\code{NULL}) sequences or removes \code{NA} values from sequences
#' in an \code{sq} object.
#'
#' @template x
#' @template by_letter
#' @template NA_letter
#' @template three-dots
#'  
#' @return An \code{\link[=sq-class]{sq}} object with the same type as the
#' input type. Sequences that do not contain any \code{NA} values are left
#' unchanged.
#' 
#' @details
#' \code{NA} may be introduced as a result of using functions like
#' \code{\link{substitute_letters}} or \code{\link{bite}}. They can also appear
#' in sequences if the user reads FASTA file using \code{\link{read_fasta}} or
#' constructs \code{sq} object from \code{\link{character}} vector with
#' \code{\link{sq}} function without \code{safe_mode} turned on - and there are
#' letters in file or strings other than specified in the alphabet.
#'
#' \code{remove_na()} is used to filter out sequences or elements that have
#' \code{NA} value(s). By default, if any letter in a sequence is \code{NA},
#' then whole sequence is replaced by empty (\code{NULL}) sequence. However, if
#' \code{by_letter} parameter is set to \code{TRUE}, then sequences are
#' only shortened by excluding \code{NA} values.
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
#'              alphabet = "dna_ext")
#'
#' # Substituting some letters with NA
#' sq_ami_sub <- substitute_letters(sq_ami, c(E = NA_character_, R = NA_character_))
#' sq_dna_sub <- substitute_letters(sq_dna, c(N = NA_character_))
#'
#' # Biting sequences out of range
#' sq_bitten <- bite(sq_ami, 1:15)
#'
#' # Printing the sequences
#' sq_ami_sub
#' sq_dna_sub
#'
#' # Removing sequences containing NA
#' remove_na(sq_ami_sub)
#' remove_na(sq_dna_sub)
#' remove_na(sq_bitten)
#'
#' # Removing only NA elements
#' remove_na(sq_ami_sub, by_letter = TRUE)
#' remove_na(sq_dna_sub, TRUE)
#' remove_na(sq_bitten, TRUE)
#'
#' @family cleaning_functions
#' @seealso \code{\link[=sq-class]{sq}}
#' @export
remove_na <- function(x, by_letter = FALSE, ...) {
  assert_flag(by_letter)
  
  UseMethod("remove_na")
}

#' @export
remove_na.default <- function(x, by_letter = FALSE, ...)
  stop("'remove_na' isn't implemented for this type of object", call. = FALSE)

#' @rdname remove_na
#' @export
remove_na.sq <- function(x, by_letter = FALSE, ...,
                         NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_remove_NA(x, by_letter, NA_letter)
}
