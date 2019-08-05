#' Read a FASTA file
#'
#' Reads a FASTA file of nucleotides or amino acids file and returns
#' a sqtibble with number of rows corresponding to the number of sequences and two
#' columns: 'name' and 'sq' giving the name of the sequence and the sequence itself.
#' @param file a \code{\link[base]{connection}} object or a \code{character} string.
#' @param type of the sequence (one of \code{ami}, \code{nuc} or \code{unt}).
#' @examples
#' \dontrun{
#' read_fasta(file = 'https://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt')
#' read_fasta("https://www.uniprot.org/uniprot/P28307.fasta")
#' }
#' @seealso \code{\link[base]{readLines}}
#' @importFrom stringi stri_detect_regex
#' @importFrom stringi stri_join
#' @export
read_fasta <- function(file, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_file_is_char(file)
  file <- .get_readable_file(file)
  
  if(.is_no_check_mode()) {
    .check_nc_is_clean_in_TRUE_FALSE(is_clean)
    .check_nc_type_in_ami_nuc(type)
    
    .nc_read_fasta(file, type, is_clean)
  } else {
    .check_is_clean_in_TRUE_FALSE_NULL(is_clean)
    .check_type_in_ami_nuc_unt_NULL(type)
    
    if (!is.null(non_standard)) {
      .nonst_read_fasta(file, type, is_clean, non_standard)
    } else {
      alph <- find_alph(file)
      .check_alph_matches_type(alph, type, is_clean)
      
      if (is.null(type)) type <- .guess_type_by_alph(alph)
      if (type != "unt" && is.null(is_clean)) {
        is_clean <- if (type == "ami") .guess_ami_is_clean(alph) else .guess_nuc_is_clean(alph) 
        .nc_read_fasta(file, type, is_clean)
      } else {
        .check_alph_length(alph)
        
        sqtibble <- read_fasta_file(file, alph)
        class(sqtibble[["sq"]]) <- c("untsq", "sq")
        attr(sqtibble[["sq"]], "alphabet") <- alph
        as_tibble(sqtibble)
      }
    }
  }
}

#' @export
write_fasta <- function(sq, name, file, nchar = 80) {
  validate_sq(sq)
  .check_name_proper_char(name)
  .check_file_is_char(file)
  .check_nchar_proper_int(nchar)
  .check_eq_lens(sq, name)
  
  sq <- .debitify_sq(sq, "char")
  char_vec <- unlist(lapply(1L:length(sq), function(i) {
    s <- sq[[i]]
    s <- lapply(split(s, floor((0:(length(s)-1))/nchar)), function(l) paste(l, collapse=""))
    paste0(">", name[i], "\n", paste(s, collapse = "\n"), "\n")
  }))
  
  writeLines(text = char_vec, con = file)
}

.nc_read_fasta <- function(file, type, is_clean) {
  sqtibble <- nc_read_fasta_file(file, type == "ami", is_clean)
  class(sqtibble[["sq"]]) <- c(if (is_clean) "clnsq" else NULL, paste0(type, "sq"), "sq")
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