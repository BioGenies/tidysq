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
#' exception of sq objects of type ami - in their alphabet there is '*' letter 
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
find_motifs <- function(sq, name, motifs) {
  UseMethod("find_motifs")
}

#' @export
find_motifs.default <- function(sq, name, motifs) {
  stop("method 'find_motifs' isn't implemented for this type of object")
}



#' @export
#' @importFrom dplyr bind_rows
#' @importFrom stringi stri_sub stri_locate_all_regex stri_count_regex
#' @importFrom tibble add_column
find_motifs.sq <- function(sq, name, motifs) {
  .validate_sq(sq)
  .check_character(name, "'name'")
  .check_eq_lens(sq, name, "'sq'", "'name'")
  .check_character(motifs, "'motifs'")
  type <- .get_sq_type(sq)
  
  sq_c <- sq
  if (type %in% c("ami", "dna", "rna"))
    motifs <- toupper(motifs)
  motifs_c <- motifs
  # Needed positive look-ahead to allow overlapping matches
  motifs_l <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs <- ifelse(motifs_l == 1,
                   motifs, 
                   paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  
  alph <- alphabet(sq)
  sq <- as.character(sq)
  
  if (type == "ami") {
    .check_motifs_proper_alph(motifs_c, "ami")
    motifs <- .replace_ami_motif(motifs)
  } else if (type == "dna") {
    .check_motifs_proper_alph(motifs_c, "dna")
    motifs <- .replace_dna_motif(motifs)
  } else if (type == "rna") {
    .check_motifs_proper_alph(motifs_c, "rna")
    motifs <- .replace_rna_motif(motifs)
  } else {
    .check_motifs_proper_alph(motifs_c, type, alph)
  }
  
  # Has to pass motifs_c so that sought column can be set
  ret_tibble <- mapply(function(motif, motif_name, motif_length) {
    ret <- stri_locate_all_regex(sq, motif, omit_no_match = TRUE)
    sequence_index <- vapply(ret, nrow, integer(1))
    ret <- add_column(
      as_tibble(do.call(rbind, ret)),
      name = rep(name, sequence_index),
      sq = rep(sq_c, sequence_index),
      sought = motif_name,
      .before = "start"
    )
    ret[["end"]] <- ret[["end"]] + rep(motif_length, nrow(ret)) - 1
    found <- stri_sub(rep(as.character(sq_c), sequence_index),
                      from = ret[["start"]],
                      to = ret[["end"]])
    ret <- add_column(
      ret,
      # TODO: replace .construct_sq_s with something... cleaner?
      found = .construct_sq_s(found, alph,
                              c(.get_sq_subclass(sq_c), if (.is_cleaned(sq_c)) "clnsq", "sq")),
      .before = "start"
    )
    # ret[!is.na(ret[, "start"]), , drop = FALSE]
  }, motifs, motifs_c, motifs_l, SIMPLIFY = FALSE)
  
  # While base::rbind would be nicer dependency-wise, it doesn't work for sq object.
  # Don't expect this to change as shown in the issue below
  # https://github.com/tidyverse/tibble/issues/34
  do.call(bind_rows, ret_tibble)
}
