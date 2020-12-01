#' Save sq to fasta file
#' 
#' Writes \code{\link[=sq-class]{sq}} objects with their names to a fasta file.
#' @param x a \code{\link[=sq-class]{sq}} object.
#' @param name a \code{\link{character}} vector of length equal to \code{sq} length.
#' @param file a \code{\link{character}} string indicating path to file to write into.
#' @param width a positive \code{\link{integer}} value informing about maximum number of 
#' characters to put in each line of file.
#' @export
write_fasta <- function(x, names, file, 
                        width = 80,
                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_class(x, "sq")
  assert_character(names, len = vec_size(x), any.missing = FALSE)
  assert_string(file)
  assert_count(width, positive = TRUE)
  assert_string(NA_letter, min.chars = 1)
  
  CPP_write_fasta(x, names, file, width, NA_letter)
}