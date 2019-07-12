#' @export
read_fasta <- function(file, type = "unt") {
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
  if (!file.exists(file)) {
    stop("'file' doesn't exists")
  }
  if (!(type %in% c("ami", "nuc", "unt"))) {
    stop("'type' needs to be one of 'ami', 'nuc', 'unt' ('unt' is default)")
  }
  
  #used from biogram
  all_lines <- readLines(file)
  s_id <- cumsum(grepl("^>", all_lines))
  all_s <- split(all_lines, s_id)
  
  s_list <- unname(sapply(all_s, function(s) paste(s[2:length(s)], collapse = "")))
  
  names_vec <- sub(">", "", sapply(all_s, function(s) s[1]), fixed = TRUE)
  
  construct_sqtibble(s_list, names_vec, type)
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
  char_vec <- unlist(lapply(1L:length(sq), function(i) {
    s <- alph[sq[[i]]]
    s <- lapply(split(s, floor((0:(length(s)-1))/nchar)), function(l) paste(l, collapse=""))
    paste0(">", name[i], "\n", paste(s, collapse = "\n"), "\n")
  }))
  
  writeLines(text = char_vec, con = file)
}