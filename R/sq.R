#' Construct sq object from character vector
#'
#' @templateVar alph_null_ok TRUE
#' 
#' @description This function allows the user to construct objects of 
#' \code{\link[=sq-class]{class sq}} from a character vector.
#' 
#' @param x [\code{character}]\cr
#'  Vector to construct \code{sq} object from.
#' @template alphabet
#' @template NA_letter
#' @template safe_mode
#' @template on_warning
#' @template ignore_case
#'
#' @return An object of \code{\link[=sq-class]{class sq}} with appropriate type.
#' 
#' @details
#' Function \code{sq} covers all possibilities of standard and non-standard
#' types and alphabets. You can check what 'type' and 'alphabet' exactly are in
#' \code{\link[=sq-class]{sq class}} documentation. There is a guide below on
#' how function operates and how the program behaves depending on arguments
#' passed and letters in the sequences.
#' 
#' \code{x} parameter should be a character vector. Each element of this vector
#' is a biological sequence. If this parameter has length 0, object of class
#' \code{sq} with 0 sequences will be created (if not specified, it will have
#' \strong{dna_bsc} type, which is a result of rules written below). If it
#' contains sequences of length 0, \code{NULL} sequences will be introduced (see
#' \emph{NULL (empty) sequences} section in \code{\link[=sq-class]{sq class}}).
#' 
#' \strong{Important note:} in all below cases word 'letter' stands for an
#' element of an alphabet. Letter might consist of more than one character, for
#' example "\code{Ala}" might be a single letter. However, if the user wants to
#' construct or read sequences with multi-character letters, one has to specify
#' all letters in \code{alphabet} parameter. Details of letters, alphabet and
#' types can be found in \code{\link[=sq-class]{sq class}} documentation.
#'
#' @section Simple guide to construct:
#' In many cases, just the \code{x} parameter needs to be specified - type of
#' sequences will be guessed according to rules described below. The user needs
#' to pay attention, however, because for short sequences type may be guessed
#' incorrectly - in this case they should specify type in \code{alphabet}
#' parameter.
#' 
#' If your sequences contain non-standard letters, where each non-standard
#' letter is one character long (that is, any character that is not an uppercase
#' letter), you also don't need to specify any parameter. Optionally, you can
#' explicitly do it by setting \code{alphabet} to \code{"unt"}.
#' 
#' In \code{\link[=safe-mode]{safe mode}} it is guaranteed that only letters
#' which are equal to \code{NA_letter} argument are interpreted as \code{NA}
#' values. Due to that, resulting alphabet might be different from the
#' \code{alphabet} argument.
#' 
#' @section Detailed guide to construct:
#' Below are listed all possibilities that can occur during the construction of
#' a \code{sq} object:
#' \itemize{
#' \item If you don't specify any other parameter than \code{x}, function will
#'  try to guess sequence type (it will check in exactly this order):
#'  \enumerate{
#'  \item If it contains only ACGT- letters, type will be set to
#'   \strong{dna_bsc}.
#'  \item If it contains only ACGU- letters, type will be set to
#'   \strong{rna_bsc}.
#'  \item If it contains any letters from 1. and 2. and additionally letters
#'   DEFHIKLMNPQRSVWY*, type will be set to \strong{ami_bsc}.
#'  \item If it contains any letters from 1. and additionally letters
#'   WSMKRYBDHVN, type will be set to \strong{dna_ext}.
#'  \item If it contains any letters from 2. and additionally letters
#'   WSMKRYBDHVN, type will be set to \strong{rna_ext}.
#'  \item If it contains any letters from previous points and additionally
#'   letters JOUXZ, type will be set to \strong{ami_ext}.
#'  \item If it contains any letters that exceed all groups mentioned above,
#'   type will be set to \strong{unt}.
#'  }
#' \item If you specify \code{alphabet} parameter as any of \code{"dna_bsc"},
#'  \code{"dna_ext"}, \code{"rna_bsc"}, \code{"rna_ext"}, \code{"ami_bsc"},
#'  \code{"ami_ext"}; then:
#'  \itemize{
#'  \item If \code{safe_mode} is \code{FALSE}, then sequences will be built
#'   with standard alphabet for given type.
#'  \item If \code{safe_mode} is \code{TRUE}, then sequences will be scanned
#'   for letters not in standard alphabet:
#'   \itemize{
#'   \item If no such letters are found, then sequences will be built with
#'    standard alphabet for given type.
#'   \item If at least one such letter is found, then sequences are built with
#'    real alphabet and with type set to \strong{unt}.
#'   }
#'  }
#' \item If you specify \code{alphabet} parameter as \code{"unt"}, then
#'  sequences are scanned for alphabet and subsequently built with obtained
#'  alphabet and type \strong{unt}.
#' \item If you specify \code{alphabet} parameter as \code{character} vector
#'  longer than 1, then type is set to \strong{atp} and alphabet is equal to
#'  letters in said parameter.
#' }
#' 
#' If \code{ignore_case} is set to \code{TRUE}, then lowercase letters are
#' turned into uppercase during their interpretation, unless type is set to
#' \strong{atp}.
#' 
#' @section Handling \strong{unt} and \strong{atp} types and \code{NA} values:
#' You can convert letters into another using \code{\link{substitute_letters}}
#' and then use \code{\link{typify}} or \code{sq_type<-} function to set type of
#' \code{sq} to \strong{dna_bsc}, \strong{dna_ext}, \strong{rna_bsc},
#' \strong{rna_ext}, \strong{ami_bsc} or \strong{ami_ext}. If your sequences
#' contain \code{NA} values, use \code{\link{remove_na}}.
#'
#' @family io_functions
#' @seealso \code{\link[=sq-class]{sq class}} \code{\link{read_fasta}}
#' \code{\link{tidysq-options}} \code{\link{substitute_letters}}
#' \code{\link{remove_na}}
#' @export
sq <- function(x,
               alphabet = NULL,
               NA_letter = getOption("tidysq_NA_letter"),
               safe_mode = getOption("tidysq_safe_mode"),
               on_warning = getOption("tidysq_on_warning"),
               ignore_case = FALSE) {
  assert_character(x, any.missing = FALSE)
  assert_flag(safe_mode)
  assert_string(NA_letter, min.chars = 1)
  assert_warning_handling(on_warning)
  assert_character(alphabet, any.missing = FALSE, min.len = 0, unique = TRUE, null.ok = TRUE)
  assert_flag(ignore_case)
  
  if (is.null(alphabet)) {
    alphabet <- obtain_alphabet(x, if (safe_mode) Inf else 4096, 
                                NA_letter, ignore_case)
    alphabet <- guess_standard_alphabet(alphabet)
  } else if (length(alphabet) == 1) {
    type <- interpret_type(alphabet)
    if (type == "unt") {
      alphabet <- obtain_alphabet(x, Inf, NA_letter, ignore_case)
    } else {
      alphabet <- get_standard_alphabet(type)
      if (safe_mode) {
        actual_alphabet <- obtain_alphabet(x, Inf, NA_letter, ignore_case)
        if (!identical(actual_alphabet, alphabet)) {
          handle_warning_message(
            "Detected letters that do not match specified type!",
            on_warning
          )
          alphabet <- actual_alphabet
        }
      }
    }
  } else {
    #TODO: safe mode should also be implemented for atp
    alphabet <- sq_alphabet(alphabet, "atp")
  }
  
  pack(x, alphabet, NA_letter, ignore_case)
}

sq_ptype <- function(str_alphabet, type)
  new_list_of(ptype = raw(0),
              alphabet = sq_alphabet(str_alphabet, type),
              class = c(type_as_class(type), "sq"))
