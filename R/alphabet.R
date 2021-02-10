#' Get alphabet of given sq object.
#' 
#' @description Returns \code{alphabet} attribute of an object.
#' 
#' @param x [\code{sq}]\cr
#'  An object to extract alphabet from.
#'  
#' @return A character vector of letters of the alphabet.
#' 
#' @details
#' Each \code{sq} object have an \strong{alphabet} associated with it. Alphabet
#' is a set of possible \strong{letters} that can appear in sequences contained
#' in object. Alphabet is kept mostly as a character vector, where each element
#' represents one \strong{letter}.
#'
#' \code{sq} objects of type \strong{ami}, \strong{dna} or \strong{rna} have
#' fixed, predefined alphabets. In other words, if two \code{sq} objects have
#' exactly the same type - \strong{ami_bsc}, \strong{dna_ext}, \strong{rna_bsc}
#' or any other combination - they are ensured to have the same alphabet.
#'
#' Below are listed alphabets for these types:
#' \itemize{
#' \item \strong{ami_bsc} - ACDEFGHIKLMNPQRSTVWY-*
#' \item \strong{ami_ext} - ABCDEFGHIJKLMNOPQRSTUVWXYZ-*
#' \item \strong{dna_bsc} - ACGT-
#' \item \strong{dna_ext} - ACGTWSMKRYBDHVN-
#' \item \strong{rna_bsc} - ACGU-
#' \item \strong{rna_ext} - ACGUWSMKRYBDHVN-
#' }
#'
#' Other types of \code{sq} objects are allowed to have different alphabets.
#' Furthermore, having an alphabet exactly identical to one of those above does
#' not automatically indicate that the type of the sequence is one of those -
#' e.g., there might be an \strong{atp} \code{sq} that has an alphabet
#' identical to \strong{ami_bsc} alphabet. To set the type, one should
#' use the \code{\link{typify}} or \code{`sq_type<-`} function.
#'
#' The purpose of co-existence of \strong{unt} and \strong{atp} alphabets is
#' the fact that although there is a standard for format of \emph{fasta} files,
#' sometimes there are other types of symbols, which do not match the standard.
#' Thanks to these types, tidysq can import files with customized alphabets.
#' Moreover, the user may want to group amino acids with similar properties
#' (e.g., for machine learning) and replace the standard alphabet with symbols
#' for whole groups. To check details, see \code{\link{read_fasta}},
#' \code{\link{sq}} and \code{\link{substitute_letters}}.
#'
#' \strong{Important note:} in \strong{atp} alphabets there is a possibility
#' of letters appearing that consist of more than one character - this
#' functionality is provided in order to handle situations like
#' post-translational modifications, (e.g., using "\code{mA}" to indicate
#' methylated alanine).
#'
#' \strong{Important note:} alphabets of \strong{atp} and \strong{unt}
#' \code{sq} objects are case sensitive. Thus, in their alphabets both
#' lowercase and uppercase characters can appear simultaneously and they are
#' treated as different letters. Alphabets of \strong{dna}, \strong{rna} and
#' \strong{ami} types are always uppercase and all functions converts other
#' parameters to uppercase when working with \strong{dna}, \strong{rna} or
#' \strong{ami} - e.g. \code{\link{\%has\%}} operator converts lower letters to
#' upper when searching for motifs in \strong{dna}, \strong{rna} or
#' \strong{ami} object.
#'
#' \strong{Important note:} maximum length of an alphabet is
#' \strong{30 letters}. The user is not allowed to read fasta files or
#' construct \code{sq} objects from character vectors that have more than 30
#' distinct characters in sequences (unless creating \strong{ami}, \strong{dna}
#' or \strong{rna} objects with \code{ignore_case} parameter set equal to
#' \code{TRUE}).
#'
#' @family alphabet_functions
#' @seealso \code{\link[=sq-class]{sq class}}
#' @export
alphabet <- function(x)
  attr(x, "alphabet")

`alphabet<-` <- function(x, value) {
  attr(x, "alphabet") <- value
  x
}

# alphabet creation ----

sq_alphabet <- function(alph, type) {
  # if exported add asserts
  new_vctr(
    alph,
    type = type,
    class = c("sq_alphabet", "character")
  )
}

sq_alphabet_ptype <- function(type)
  sq_alphabet(character(), type)
#' Get standard alphabet for given type.
#'
#' @description Returns \code{alphabet} attribute of an object.
#'
#' @param type [\code{character(1)}]\cr
#'  The name of standard sq type - one of \code{"dna_bsc"}, \code{"dna_ext"},
#'  \code{"rna_bsc"}, \code{"rna_ext"}, \code{"ami_bsc"} and \code{"ami_ext"}.
#'
#' @return An \code{sq_alphabet} object related to passed sq type.
#'
#' @details
#' Each of standard sq types has exactly one predefined alphabet. It allows
#' \pkg{tidysq} to package to optimize type-specific operations like
#' \code{\link{complement}()} or \code{\link{translate}()}. This function
#' enables the user to access \code{alphabet} attribute common for all \code{sq}
#' objects of given type.
#'
#' For list of letters specific to any of these standard alphabets, see
#' \code{\link{alphabet}()}.
#'
#' @family alphabet_functions
#' @export
get_standard_alphabet <- function(type) {
  type <- interpret_type(type)
  CPP_get_standard_alphabet(type)
}

# Scans ProtoSq object to determine alphabet (understood as unique encountered letters).
obtain_alphabet <- function(x, sample_size = 4096, 
                            NA_letter = getOption("tidysq_NA_letter"),
                            ignore_case = FALSE) {
  CPP_obtain_alphabet(x, sample_size, NA_letter, ignore_case)
}

# Finds smallest standard alphabet that contains all letters from alph.
guess_standard_alphabet <- function(alph,
                                    NA_letter = getOption("tidysq_NA_letter")) {
  CPP_guess_standard_alph(alph, NA_letter)
}

# utility methods ----

`[.sq_alphabet` <- function(x, i,
                            NA_letter = getOption("tidysq_NA_letter")) {
  ret <- vec_data(x)[i]
  ret[i == (2 ^ size(x) - 1)] <- NA_letter
  ret
}

# Returns number of bits that each letter uses.
size <- function(alph) {
  ceiling(log2(length(alph) + 1))
}
