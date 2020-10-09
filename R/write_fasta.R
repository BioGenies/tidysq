#' Save sq to fasta file
#' 
#' Writes \code{\link{sq}} objects with their names to a fasta file.
#' @param sq a \code{\link{sq}} object.
#' @param name a \code{\link{character}} vector of length equal to \code{sq} length.
#' @param file a \code{\link{character}} string indicating path to file to write into.
#' @param nchar a positive \code{\link{integer}} value informing about maximum number of 
#' characters to put in each line of file.
#' @export
write_fasta <- function(sq, name, file, nchar = 80) {
  assert_class(sq, "sq")
  assert_character(name, len = vec_size(sq), any.missing = FALSE)
  assert_string(file)
  assert_count(nchar, positive = TRUE)
  
  sq <- .unpack_from_sq(sq, "char")
  char_vec <- unlist(lapply(1L:length(sq), function(i) {
    s <- sq[[i]]
    s <- lapply(split(s, floor((0:(length(s) - 1))/nchar)), function(l) paste(l, collapse = ""))
    paste0(">", name[i], "\n", paste(s, collapse = "\n"), "\n")
  }))
  
  writeLines(text = char_vec, con = file)
}