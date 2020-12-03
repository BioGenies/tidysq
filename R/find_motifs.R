#' Find given motifs
#'
#' @templateVar name_null_ok FALSE
#' 
#' @description Find all given motifs in all sequences and return their 
#' positions.
#' 
#' @template x
#' @template name
#' @param motifs [\code{character}]\cr
#'  Motifs to be searched for.
#' @template NA_letter
#' @template three-dots
#' 
#' @return A \code{\link[tibble]{tibble}} with following columns:
#'  \item{name}{name of the sequence}
#'  \item{sq}{sequence}
#'  \item{sought}{sought motif}
#'  \item{found}{motif found in a sequence, may differ from sought if a motif
#'  contained ambiguous letters}
#'  \item{start}{position of motif start}
#'  \item{end}{position of motif end}
#' 
#' @details This function allows search of a given motif or motifs in the \code{sq} 
#' object. It returns all motifs found with their start and end positions 
#' within a sequence.
#' 
#' @section Allowed and forbidden letters and characters details:
#' Note if a sq object contains characters: ^$?=()\.|+*{}[] in its alphabet, 
#' search for motifs cannot be performed and an error will be displayed (with 
#' exception of sq objects of type ami - there is '*' letter in their alphabet
#' and it can be contained in sought motif). To search for motifs with those 
#' characters, you have to replace them first using 
#' \code{\link{substitute_letters}}. 
#' 
#' If sq objects of type \strong{ami}, \strong{dna} and \strong{rna}, motifs have to
#' consist of upper case letters from amino acid, DNA and RNA alphabets respectively.
#' Use of lower case letters will return an error. Two additional characters 
#' are allowed: '^' and '$' indicating the beginning and the end of a sequence 
#' respectively. Moreover, notice that '*' character may be used in amino acid 
#' motifs, as it is a part of the amino acid alphabet. If a motif contains 
#' ambiguous letters, all possible matches will be searched for. For example the 
#' amino acid motif "MAJ" (where "J" is an ambiguous letter indicating L or I) will 
#' find the motifs: "MAJ", "MAL" and "MAI". 
#' 
#' Detailed list of all letters corresponding to each ambiguous letter may be found at
#' \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}}.
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link{substitute_letters}} \code{\link{\%has\%}}
#' 
#' @export
find_motifs <- function(x, name, motifs, ...) {
  assert_character(name, len = vec_size(x))
  assert_character(motifs, any.missing = FALSE)
  
  UseMethod("find_motifs")
}

#' @export
find_motifs.default <- function(x, name, motifs, ...)
  stop("method 'find_motifs' isn't implemented for this type of object")

#' @rdname find_motifs
#' @export
#' @importFrom stringi stri_sub
#' @importFrom tibble as_tibble
find_motifs.sq <- function(x, name, motifs, ...,
                           NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  assert_alph_regex_friendly(alphabet(x))
  
  ret <- CPP_find_motifs(x, name, motifs, NA_letter)
  as_tibble(ret)
}
