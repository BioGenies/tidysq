#' Read a FASTA file
#'
#' Reads a FASTA file of nucleotides or amino acids file and returns
#' a \code{\link[tibble]{tibble}} with number of rows corresponding to the number of sequences 
#' and two columns: 'name' and 'sq' giving the name of the sequence and the sequence itself.
#' @param file a \code{\link{character}} string indicating path to file or url.
#' @inheritParams construct_sq
#' @return 
#' A \code{\link[tibble]{tibble}} with number of rows corresponding to the number of sequences 
#' and two columns: 'name' and 'sq' giving the name of the sequence and the sequence itself.
#' @details 
#' All rules of creating sq objects are the same as in \code{\link{construct_sq}}.
#' 
#' Functions \code{read_fasta_ami} and \code{read_fasta_nuc} are wrappers around 
#' \code{read_fasta} with specified \code{type} parameter - accordingly "ami" or "nuc". You
#' can also pass "is_clean" parameter to those functions, but you cannot pass "non_standard".
#' 
#' @examples
#' fasta_file <- system.file(package = "tidysq", 
#'                      "sample_fasta/sample_ami.fasta")
#' read_fasta(fasta_file)
#' read_fasta_ami(fasta_file)
#' \dontrun{
#' read_fasta("https://www.uniprot.org/uniprot/P28307.fasta")
#' }
#' @seealso \code{\link[base]{readLines}} \code{\link{construct_sq}}
#' @importFrom stringi stri_detect_regex stri_join
#' @export
read_fasta <- function(file, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_character(file, "'file'", single_elem = TRUE)
  file <- .get_readable_file(file)
  
  if (.is_fast_mode()) {
    .check_logical(is_clean, "'is_clean'", single_elem = TRUE)
    .check_type(type)
    .nc_read_fasta(file, type, is_clean)
  } else {
    .check_logical(is_clean, "'is_clean'", single_elem = TRUE, allow_null = TRUE)
    .check_type(type, allow_unt = TRUE, allow_null = TRUE)
    
    if (!is.null(non_standard)) {
      .nonst_read_fasta(file, type, is_clean, non_standard)
    } else {
      alph <- find_alph(file)
      if (!is.null(type) && type %in% c("ami", "dna", "rna")) alph <- toupper(alph)
      .check_alph_matches_type(alph, type, is_clean)
      
      if (is.null(type)) {
        type_clean <- .guess_type_subtype_by_alph(alph)
        type <- type_clean[["type"]]
        if (is.null(is_clean) && type != "unt") is_clean <- type_clean[["is_clean"]]
      } else if (type != "unt" && is.null(is_clean)) {
        if      (type == "ami") is_clean <- .guess_ami_is_clean(alph)
        else if (type == "dna") is_clean <- .guess_dna_is_clean(alph)
        else if (type == "rna") is_clean <- .guess_rna_is_clean(alph)
      } 
      if (type != "unt") {
        .nc_read_fasta(file, type, is_clean)
      } else {
        .check_alph_length(alph)
        
        sqtibble <- read_fasta_file(file, alph)
        class(sqtibble[["sq"]]) <- c("untsq", "sq", "list")
        attr(sqtibble[["sq"]], "alphabet") <- alph
        as_tibble(sqtibble)
      }
    }
  }
}

#' @rdname read_fasta
#' @export
read_fasta_ami <- function(file, is_clean = NULL) {
  read_fasta(file, type = "ami", is_clean)
}

#' @rdname read_fasta
#' @export
read_fasta_dna <- function(file, is_clean = NULL) {
  read_fasta(file, type = "dna", is_clean)
}

#' @rdname read_fasta
#' @export
read_fasta_rna <- function(file, is_clean = NULL) {
  read_fasta(file, type = "rna", is_clean)
}

.nc_read_fasta <- function(file, type, is_clean) {
  sqtibble <- nc_read_fasta_file(file, "ami", is_clean)
  class(sqtibble[["sq"]]) <- c(if (is_clean) "clnsq" else NULL, paste0(type, "sq"), "sq", "list")
  attr(sqtibble[["sq"]], "alphabet") <- .get_standard_alph(type, is_clean)
  as_tibble(sqtibble)
}

.nonst_read_fasta <- function(file, type, is_clean, non_standard) {
  all_lines <- readLines(file)
  s_id <- cumsum(stri_detect_regex(all_lines, "^>"))
  all_s <- split(all_lines, s_id)
  s_list <- unname(sapply(all_s, function(s) stri_join(s[-1], collapse = "")))
  sq <- construct_sq(s_list, type, is_clean, non_standard)
  names_vec <- stri_sub(sapply(all_s, function(s) s[1]), 2)
  
  tibble(name = names_vec, sq = sq)
}