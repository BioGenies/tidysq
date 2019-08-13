#' Read a FASTA file
#'
#' Reads a FASTA file of nucleotides or amino acids file and returns
#' a sqtibble with number of rows corresponding to the number of sequences and two
#' columns: 'name' and 'sq' giving the name of the sequence and the sequence itself.
#' @param file a \code{\link[base]{connection}} object or a \code{character} string.
#' @param type of the sequence (one of \code{ami}, \code{nuc} or \code{unt}).
#' @param is_clean \code{logical}, if \code{TRUE}, it is assumed that all 
#' sequences do not contain ambiguous letters. 
#' @param non_standard \code{character} vector of non-standard letters/
#' @examples
#' read_fasta(system.file(package = "tidysq", 
#'                      "sample_fasta/sample_ami.fasta"))
#' \dontrun{
#' read_fasta("https://www.uniprot.org/uniprot/P28307.fasta")
#' }
#' @seealso \code{\link[base]{readLines}}
#' @importFrom stringi stri_detect_regex
#' @importFrom stringi stri_join
#' @export
read_fasta <- function(file, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_character(file, "'file'", single_elem = TRUE)
  file <- .get_readable_file(file)
  
  if (.is_no_check_mode()) {
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
      if (!is.null(type) && type %in% c("ami", "nuc")) alph <- toupper(alph)
      .check_alph_matches_type(alph, type, is_clean)
      
      if (is.null(type)) {
        type_clean <- .guess_type_subtype_by_alph(alph)
        type <- type_clean[["type"]]
        if (is.null(is_clean) && type != "unt") is_clean <- type_clean[["is_clean"]]
      } else if (type != "unt" && is.null(is_clean)) {
        is_clean <- if (type == "ami") .guess_ami_is_clean(alph) else .guess_nuc_is_clean(alph) 
        
      } 
      if (type != "unt") {
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
  .check_character(name, "'name'")
  .check_character(file, "'file'", single_elem = TRUE)
  .check_integer(nchar, "'nchar'", single_elem = TRUE, allow_negative = FALSE, allow_zero = FALSE)
  .check_eq_lens(sq, name, "'sq'", "'name'")
  
  sq <- .debitify_sq(sq, "char")
  char_vec <- unlist(lapply(1L:length(sq), function(i) {
    s <- sq[[i]]
    s <- lapply(split(s, floor((0:(length(s) - 1))/nchar)), function(l) paste(l, collapse = ""))
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