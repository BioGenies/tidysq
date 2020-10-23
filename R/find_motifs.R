#' Find given motifs
#' 
#' @description Find all given motifs in all sequences and return their 
#' positions.
#' 
#' @inheritParams reverse
#' @param name a non-\code{NULL} \code{character} vector without \code{\link{NA}} values, 
#' containing names of the sequences in the sq. It has to be of the same length 
#' as the \code{sq}. 
#' @param motifs a \code{character} vector of motifs to be searched for.
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
#' @examples 
#' # Creating objects to work on:
#' sq_ami <- construct_sq(c("AGNTYIKFGGAYTI", "MATEGILIAADGYTWIL", 
#'                          "MIPADHICAANGIENAGIK"), type = 'ami')
#' sqtbl <- read_fasta(system.file(package = "tidysq", "example_aa.fasta"), 
#'                     type = "ami")
#' 
#' # Find motif of two alanines followed by aspartic acid or asparagine 
#' # ('AAB' motif will match 'AAB', 'AAD' and 'AAN'):
#' find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "AAB")
#' 
#' # Find motif 'VHH' at the beginning of sequences:
#' find_motifs(sqtbl[["sq"]], sqtbl[["name"]], "^VHH")
#' 
#' # Find motif 'DPGS' at the end of sequences:
#' find_motifs(sqtbl[["sq"]], sqtbl[["name"]], "DPGS$")
#' 
#' # Find multiple motifs:
#' find_motifs(sqtbl[["sq"]], sqtbl[["name"]], c("^LIV", "XXKK", "EN$"))
#' 
#' @seealso \code{\link{sq}} \code{\link{substitute_letters}} \code{\link{\%has\%}}
#' 
#' @export
find_motifs <- function(x, name, motifs) {
  assert_character(name, len = vec_size(x))
  assert_character(motifs, any.missing = FALSE)
  
  UseMethod("find_motifs")
}

#' @export
find_motifs.default <- function(x, name, motifs)
  stop("method 'find_motifs' isn't implemented for this type of object")

#' @export
#' @importFrom stringi stri_sub
find_motifs.sq <- function(x, name, motifs) {
  assert_alph_regex_friendly(alphabet(x))
  
  motif_lengths <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs_regex <- ifelse(motif_lengths == 1,
                         motifs,
                         paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  
  .find_motifs_sq(x, name, motifs, motifs_regex, motif_lengths)
}

#' @export
#' @importFrom stringi stri_sub
find_motifs.dnasq <- function(x, name, motifs) {
  motifs <- toupper(motifs)
  assert_motifs_for_type(motifs, "dna")
  
  motif_lengths <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs_regex <- ifelse(motif_lengths == 1,
                         motifs,
                         paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  motifs_regex <- .replace_dna_motif(motifs_regex)
  
  .find_motifs_sq(x, name, motifs, motifs_regex, motif_lengths)
}

#' @export
#' @importFrom stringi stri_sub
find_motifs.rnasq <- function(x, name, motifs) {
  motifs <- toupper(motifs)
  assert_motifs_for_type(motifs, "rna")
  
  motif_lengths <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs_regex <- ifelse(motif_lengths == 1,
                         motifs,
                         paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  motifs_regex <- .replace_rna_motif(motifs_regex)
  
  .find_motifs_sq(x, name, motifs, motifs_regex, motif_lengths)
}

#' @export
#' @importFrom stringi stri_sub
find_motifs.amisq <- function(x, name, motifs) {
  motifs <- toupper(motifs)
  assert_motifs_for_type(motifs, "ami")
  
  motif_lengths <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs_regex <- ifelse(motif_lengths == 1,
                         motifs,
                         paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  motifs_regex <- .replace_ami_motif(motifs_regex)
  
  .find_motifs_sq(x, name, motifs, motifs_regex, motif_lengths)
}

#' @importFrom dplyr bind_rows
#' @importFrom stringi stri_sub stri_locate_all_regex stri_count_regex
#' @importFrom tibble add_column
.find_motifs_sq <- function(x, name, motifs, motifs_regex, motif_lengths) {
  sq_character <- as.character(x)
  
  # Has to pass motifs_c so that sought column can be set
  ret_tibble <- mapply(function(motif_name, motif_regex, motif_length) {
    ret <- stri_locate_all_regex(sq_character, motif_regex, omit_no_match = TRUE)
    sequence_index <- vapply(ret, nrow, integer(1))
    ret <- add_column(
      as_tibble(do.call(rbind, ret)),
      name = rep(name, sequence_index),
      sq = rep(x, sequence_index),
      sought = motif_name,
      .before = "start"
    )
    ret[["end"]] <- ret[["end"]] + rep(motif_length, nrow(ret)) - 1
    found <- stri_sub(rep(as.character(x), sequence_index),
                      from = ret[["start"]],
                      to = ret[["end"]])
    ret <- add_column(
      ret,
      found = sq(found, get_sq_type(x)),
      .before = "start"
    )
    # ret[!is.na(ret[, "start"]), , drop = FALSE]
  }, motifs, motifs_regex, motif_lengths, SIMPLIFY = FALSE)
  
  # While base::rbind would be nicer dependency-wise, it doesn't work for sq object.
  # Don't expect this to change as shown in the issue below
  # https://github.com/tidyverse/tibble/issues/34
  do.call(bind_rows, ret_tibble)
}
