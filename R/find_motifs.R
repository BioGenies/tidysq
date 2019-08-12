#' Find given motifs
#' 
#' @description Find all given motifs in all sequences and return their 
#' positions 
#' 
#' @param sq \code{\link{sq}} object to be tested
#' @param name non-NULL \code{character} vector without NA's, containing
#' names of the sequences in the sq. It have to be of the same length 
#' as the sq. 
#' @param motifs \code{character} vector of motifs to be searched for
#' 
#' @return a tibble with number of rows the same as the length of sq and
#' following columns:
#'  \item{name}{name of the sequence}
#'  \item{sq}{sequence}
#'  \item{sought}{sought motif}
#'  \item{found}{motif found in a sequence, may differ from sought if a motif
#'  contained ambiguous letters}
#'  \item{start}{position of motif start}
#'  \item{end}{position of motif end}
#' 
#' @details This function allows search of a given motif or motifs in the sq 
#' object. It returns all found motifs with their start and end positions 
#' within a sequence.
#' 
#' Note if a sq object contains characters: ^$?=()\.|+*{}[] in its alphabet, 
#' search for motifs cannot be performed and an error will be displayed (with 
#' exception of sq objects of type ami - in their alphabet there is "*" letter 
#' and it can be contained in sought motif). To search for motifs with those 
#' characters, you have to replace them first using 
#' \code{\link{substitute_letters}}. 
#' 
#' In case of sq objects of type 'ami' and 'nuc', motifs have to consist of 
#' upper case letters from amino acid and nucleotide alphabet respectively. 
#' Use of lower case letters will return an error. Two additional characters 
#' are allowed: '^' and '$' indicating the beginning and the end of a sequence 
#' respectively. Moreover, notice that '*' character may be used in amino acid 
#' motifs, as it is a part of the amino acid alphabet. If a motif contains 
#' ambiguous letters, all possible matches will be searched for, e.g., amino 
#' acid motif "MAJ" (where "J" is an ambiguous letter indicating L or I) will 
#' find motifs: "MAJ", "MAL" and "MAI". 
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
#' @seealso sq substitute_letters \%has\%
#' 
#' @importFrom stringi stri_sub stri_locate_all_regex stri_count_regex
#' @export
find_motifs <- function(sq, name, motifs) {
  validate_sq(sq)
  .check_character(name, "'name'")
  .check_eq_lens(sq, name, "'sq'", "'name'")
  .check_character(motifs, "'motifs'")
  type <- .get_sq_type(sq)
  
  sq_c <- sq
  motifs_c <- motifs
  motifs_l <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs <- ifelse(motifs_l == 1,
                   motifs, 
                   paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  motifs <- strsplit(motifs, "")
  
  alph <- .get_alph(sq)
  
  if (type == "ami") {
    motifs <- lapply(motifs, toupper)
    .check_motifs_proper_alph(motifs_c, "ami")
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "*", "\\*"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "B", "[BDN]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "J", "[JIL]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "X", "[A-Z]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "Z", "[ZEQ]"))
    motifs <- sapply(motifs, function(motif) paste(motif, collapse = ""))
  } else if (type == "nuc") {
    motifs <- lapply(motifs, toupper)
    .check_motifs_proper_alph(motifs_c, "nuc")
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "W", "[WATU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "S", "[SCG]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "M", "[MAC]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "K", "[KGTU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "R", "[RAG]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "Y", "[YCTU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "B", "[BCTGU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "D", "[DATGU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "H", "[HACTU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "V", "[VACG]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "N", "[ACTGUWSMKRYBDHVN]"))
    motifs <- sapply(motifs, function(motif) paste(motif, collapse = ""))
    
  } else .check_motifs_proper_alph(motifs_c, type, alph)

  sq <- as.character(sq)
  n <- length(sq)
  m <- length(motifs)
  match_inds <- lapply(1:m, function(i) {
    ret <- stri_locate_all_regex(sq, motifs[i])
    ret <- cbind(do.call(rbind, ret), s_ind = unlist(lapply(1:n, function(i) rep(i, nrow(ret[[i]])))))
    ret <- ret[!is.na(ret[,"start"]), , drop = FALSE]
    ret[, "end"] <- ret[, "end"] + rep(motifs_l[i], nrow(ret)) - 1
    ret
  })
  
  matched_inds <- cbind(do.call(rbind, match_inds)) 
  sought <-  unlist(lapply(1:m, function(j) rep(motifs_c[j], nrow(match_inds[[j]]))))
  
  sq_col <- sq_c[matched_inds[, "s_ind"]]
  nm_col <- name[matched_inds[, "s_ind"]]
  found <- stri_sub(sq[matched_inds[, "s_ind"]], from = matched_inds[, "start"], to = matched_inds[, "end"])
  sq_col <- .set_class_alph(sq_col, sq_c)
  tibble(name = nm_col, 
         sq = sq_col, 
         sought = sought, 
         found = .construct_sq_s(found, alph, class(sq_c)), 
         start = matched_inds[, "start"], 
         end = matched_inds[, "end"])
}