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
#' @export
read_fasta <- function(file, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_file_is_char(file)
  if(.is_no_check_mode()) {
    .check_nc_is_clean_in_TRUE_FALSE(is_clean)
    .check_nc_type_in_ami_nuc(type)
    sqtibble <- nc_read_fasta_file(normalizePath(file), type == "ami", is_clean)
    class(sqtibble[["sq"]]) <- c(if (is_clean) "clnsq" else NULL, paste0(type, "sq"), "sq")
    attr(sqtibble[["sq"]], "alphabet") <- .get_standard_alph(type, is_clean)
    as_tibble(sqtibble)
  } else {
    .check_is_clean_in_TRUE_FALSE_NULL(is_clean)
    
    #used from biogram
    all_lines <- readLines(file)
    s_id <- cumsum(grepl("^>", all_lines))
    all_s <- split(all_lines, s_id)
    
    s_list <- unname(sapply(all_s, function(s) paste(s[2:length(s)], collapse = "")))
    sq <- construct_sq(s_list, type, is_clean, non_standard)
    
    names_vec <- sub(">", "", sapply(all_s, function(s) s[1]), fixed = TRUE)
    
    tibble(name = names_vec, sq = sq)
  }
}

#' @export
write_fasta <- function(sq, name, file, nchar = 80) {
  validate_sq(sq)
  if (missing(name) ||
      !is.character(name) ||
      any(is.null(name) | is.na(name)) ||
      (length(name) != length(sq))) {
    stop("'name' has to be a non-NULL character vector without NA's, of lenght equal to length of sq")
  }
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
  if (!is.numeric(nchar) ||
      (floor(nchar) != nchar) ||
      (length(nchar) != 1) ||
      is.na(nchar) || 
      is.nan(nchar)|| 
      !is.finite(nchar) ||
      nchar <= 0) {
    stop("'nchar' has to be positive integer indicating max number of characters in single line in file")
  }
  
  alph <- .get_alph(sq)
  sq <- .debitify_sq(sq, alph)
  char_vec <- unlist(lapply(1L:length(sq), function(i) {
    s <- sq[[i]]
    s <- lapply(split(s, floor((0:(length(s)-1))/nchar)), function(l) paste(l, collapse=""))
    paste0(">", name[i], "\n", paste(s, collapse = "\n"), "\n")
  }))
  
  writeLines(text = char_vec, con = file)
}

read_fasta_nc <- function(file, type, is_clean = TRUE) {
  if (missing(type) ||
      !type %in% c("nuc", "ami")) {
    stop("in no_check mode 'type' needs to be one of 'nuc', 'ami'")
  }
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
  if (!file.exists(file)) {
    stop("'file' doesn't exists")
  }
  if (!is_clean %in% c(TRUE, FALSE)) {
    stop("'is_clean' has to be TRUE or FALSE")
  }
  
  #used from biogram
  all_lines <- readLines(file)
  s_id <- cumsum(grepl("^>", all_lines))
  all_s <- split(all_lines, s_id)
  
  s_list <- unname(sapply(all_s, function(s) paste(s[2:length(s)], collapse = "")))
  sq <- construct_sq_nc(s_list, type, is_clean)
  
  names_vec <- sub(">", "", sapply(all_s, function(s) s[1]), fixed = TRUE)
  
  construct_sqtibble(sq, names_vec, type)
}