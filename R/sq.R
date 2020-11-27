#' Construct sq object from character vector
#' 
#' @description This function allows the user to construct objects of 
#' \code{\link[=sq-class]{class sq}} from a character vector.
#' 
#' @param x [\code{character}]\cr
#'  Vector to construct object from.
#' @param alphabet [\code{character}]\cr
#'  If provided value is a single string, it will be interpreted as type (see
#'  details). If provided value has length greater than one, it will be treated
#'  as atypical alphabet for \code{sq} obejct and \code{sq} type will be
#'  \code{atp}. If provided value is \code{NULL}, type guessing will be
#'  performed (see details).
#' @template NA_letter
#' @template safe_mode
#' @template on_warning
#' @template ignore_case
#'
#' @return object of \code{\link[=sq-class]{class sq}} with appropriate type (one of: \strong{ami},
#' \strong{dna}, \strong{rna}, \strong{unt}, \strong{atp}).
#' 
#' @details
#' Function \code{construct_sq} covers all possibilities of standard and non-standard types and 
#' alphabets. You can check what 'type' and 'alphabet' exactly are in \code{\link[=sq-class]{sq class}}
#' documentation. Below there is a guide how function operates and how the program behaves 
#' according to the given arguments and the letters in the sequences.
#' 
#' Functions \code{construct_sq_ami}, \code{construct_sq_dna} and \code{construct_sq_rna} are
#' wrappers around \code{construct_sq} with specified \code{type} parameter - accordingly
#' "ami", "dna" or "rna". You can also pass "is_clean" parameter to those functions, but you
#' cannot pass "non_standard".
#' 
#' \code{sq} parameter should be a character vector. Each element of this vector is a biological 
#' sequence. If this parameter has length 0, object of class \code{sq} with 0 sequences will be 
#' created (if not specified, it will have \strong{dna} \strong{cln} type, which is a result of 
#' rules written below). If it contains sequences of length 0, \code{\link[=sq-class]{NULL}} sequences
#' will be introduced (see \emph{NULL (empty) sequences} section in \code{\link[=sq-class]{sq class}}).
#' 
#' \strong{Important note:} in all below cases word 'letter' stands for an element of an alphabet.
#' Letter might consist of more than one character, for example "Ala" might be a single letter.
#' However, if you want to construct or read sequences with multi-character letters, one has 
#' to specify \code{non_standard} parameter. Details of letters, alphabet and types can be 
#' found in \code{\link[=sq-class]{sq class}} documentation.
#'
#' @section Simple guide to construct :
#' In most cases, just the \code{sq} parameter needs to be specified - type of sequences
#' will be guessed accordingly to rules described below. You need to pay attention, however, 
#' because for short sequences type may be guessed incorrectly - in this case you should
#' specify \code{type} and/or \code{is_clean}.
#' 
#' If your sequences contain non-standard letters, where each non-standard letter is one
#' character long, you also don't need to specify any parameter. Optionally, you can explicitly
#' do it by setting \code{type} to "unt".
#' 
#' If you want to construct sequences with multicharacter letters, you have to specify 
#' \code{non_standard} parameter, where you have to provide all non-standard letters longer
#' than one character.
#' 
#' In \code{\link[=fast-mode]{fast mode}} you have to specify both \code{type} and \code{is_clean} 
#' parameters. You cannot specify \code{non_standard} parameter in this mode. All letters
#' outside specified alphabet will be red as \code{NA} values.
#' 
#' @section Detailed guide to construct :
#' Below all possibilities that can occur during the construction of a \code{sq} object are listed.
#' 
#' In normal mode (no \code{\link[=fast-mode]{fast mode}}):
#' \itemize{
#' \item If you don't specify any other parameter than \code{sq}, function will try to guess
#' sequence type and if it is clean (it will check in exactly this order):
#' \enumerate{
#' \item If it contains only ACGT- letters, either lowercase or uppercase, type will be set
#' to \strong{dna} \strong{cln}.
#' \item If it contains only ACGU- letters, either lowercase or uppercase, type will be set
#' to \strong{dna} \strong{cln}.
#' \item If it contains any letters from 1. and 2. and additionally letters DEFHIKLMNPQRSVWY*,
#' either lowercase or uppercase, type will be set to \strong{ami} \strong{cln}.
#' \item If it contains any letters from 1. and additionally letters WSMKRYBDHVN, either
#' lowercase or uppercase, and does not contain any other letter, type will be set to 
#' \strong{dna} without \strong{cln} subtype.
#' \item If it contains any letters from 2. and additionally letters WSMKRYBDHVN, either
#' lowercase or uppercase, and does not contain any other letter, type will be set to 
#' \strong{rna} without \strong{cln} subtype.
#' \item If it contains any letters from previous points and additionally letters JOUXZ, type will
#' be set to \strong{ami} without \strong{cln} subtype.
#' \item If it contains any letters that exceed all groups mentioned above, type will be set
#' to "unt".
#' }
#' \item If you specify \code{type} parameter as "ami", "dna" or "rna" (and do not specify neither 
#' \code{is_clean} nor \code{non_standard}) type will be checked during construction: 
#' \enumerate{
#' \item If all letters in sequences fit the clean alphabet of given type, type will be set to 
#' given type with \strong{cln} subtype. 
#' \item If all letters in sequences fit the unclean alphabet of given type, type will be set to 
#' given type without \strong{cln} subtype.
#' \item If at least one sequence contains at least one letter that is not an element of an  
#' unclean alphabet of the provided type, an error will be thrown. 
#' }
#' \item If you specify both \code{type} and \code{is_clean}, function checks if letters
#' in sequences matches exactly specified alphabet (with capitalisation accuracy). If they do, 
#' type will be set to it. Otherwise, an error will be thrown.
#' \item If you specify \code{type} as "unt" and neither \code{is_clean} nor \code{non_standard},
#' type will be set to \strong{unt}. Letters won't be converted to uppercase, alphabet will
#' consist of all letters found in sequences. 
#' \item If you do not specify neither \code{type} nor \code{is_clean} and specify 
#' \code{non_standard} parameter, which should be character vector where each element is at least 
#' two characters long, all strings as specified will be detected in sequences and treated as 
#' letters in constructed \strong{atp} \code{sq}.
#' \item All other combinations of parameters are incorrect.
#' }
#' 
#' In \code{\link[=fast-mode]{fast mode}} you have to specify \code{type} (it has to have one of
#' "ami", "dna" or "rna" values) and \code{is_clean} (\code{TRUE} or \code{FALSE}). You cannot
#' specify \code{non_standard}. All letters that aren't elements of destination alphabet (with
#' a letter size accuracy) will be treated as \code{NA} values.
#' 
#' @section Handling atp and unt sequences and NA values:
#' You can convert letters into another using \code{\link{substitute_letters}} and then you
#' can use \code{\link{typify}} function to set of \code{sq} to \strong{ami}, \strong{dna}
#' or \strong{rna}. If your sequences contain \code{NA} values, use \code{\link{remove_na}}
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link{read_fasta}} \code{\link{tidysq-options}}
#' \code{\link{fast-mode}} \code{\link{substitute_letters}} \code{\link{remove_na}}
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
