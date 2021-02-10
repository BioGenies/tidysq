#' Find given motifs
#'
#' @templateVar name_null_ok FALSE
#' 
#' @description Finds all given motifs in sequences and returns their positions.
#' 
#' @template x
#' @template name
#' @param motifs [\code{character}]\cr
#'  Motifs to be searched for.
#' @template NA_letter
#' @template three-dots
#' 
#' @return A \code{\link[tibble]{tibble}} with following columns:
#'  \item{name}{name of the sequence in which a motif was found}
#'  \item{sought}{sought motif}
#'  \item{found}{found subsequence, may differ from sought if the motif
#'   contained ambiguous letters}
#'  \item{start}{position of first element of found motif}
#'  \item{end}{position of last element of found motif}
#' 
#' @details
#' This function allows search of a given motif or motifs in the \code{sq}
#' object. It returns all motifs found with their start and end positions within
#' a sequence.
#' 
#' @template motif_details
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("ATGCAGGA", "GACCGNBAACGAN", "TGACGAGCTTAG"),
#'              alphabet = "dna_bsc")
#' sq_ami <- sq(c("AGNTYIKFGGAYTI", "MATEGILIAADGYTWIL", "MIPADHICAANGIENAGIK"),
#'              alphabet = "ami_bsc")
#' sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
#'              alphabet = c("mA", "mY", "nbA", "nsA"))
#' sq_names <- c("sq1", "sq2", "sq3")
#'
#' # Finding motif of two alanines followed by aspartic acid or asparagine
#' # ("AAB" motif matches "AAB", "AAD" and "AAN"):
#' find_motifs(sq_ami, sq_names, "AAB")
#'
#' # Finding "C" at fourth position:
#' find_motifs(sq_dna, sq_names, "^NNNC")
#'
#' # Finding motif "I" at second-to-last position:
#' find_motifs(sq_ami, sq_names, "IX$")
#'
#' # Finding multiple motifs:
#' find_motifs(sq_dna, sq_names, c("^ABN", "ANCBY", "BAN$"))
#'
#' # Finding multicharacter motifs:
#' find_motifs(sq_atp, sq_names, c("nsA", "mYmY$"))
#'
#' @family bio_functions
#' @export
find_motifs <- function(x, name, motifs, ...) {
  assert_character(name, len = vec_size(x), unique = TRUE)
  assert_character(motifs, any.missing = FALSE)
  
  UseMethod("find_motifs")
}

#' @export
find_motifs.default <- function(x, name, motifs, ...)
  stop("method 'find_motifs' isn't implemented for this type of object")

#' @rdname find_motifs
#' @export
#' @importFrom tibble as_tibble
find_motifs.sq <- function(x, name, motifs, ...,
                           NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  assert_alph_no_special_chars(alphabet(x))
  
  ret <- CPP_find_motifs(x, name, motifs, NA_letter)
  as_tibble(ret)
}
