#' Apply function to each sequence
#' 
#' Applies given function to each sequence. Sequences are passed to function as character vectors
#' (or numeric, if type of \code{sq} is \strong{enc}) or single character strings, depending on 
#' parameter.
#' 
#' @template x
#' @param fun [\code{function(1)}]\cr
#'  A function to apply to each sequence in \code{sq} object; it should
#'  take a character vector, numeric vector or single character string as an input.
#' @template three-dots
#' @param single_string [\code{logical(1)}]\cr
#'  A value indicating in which form sequences should be
#'  passed to the function \code{fun}; if \code{FALSE} (default), they will be treated as character
#'  vectors, if \code{TRUE}, they will be pasted into a single string.
#' @template NA_letter
#' 
#' @return A list of values returned by function for each sequence in corresponding order.
#' 
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("ATGCAGGA", "GACCGNBAACGAN", "TGACGAGCTTA"),
#'              alphabet = "dna_bsc")
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_unt <- sq(c("ATGCAGGA?", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#'
#' # Counting how may "A" elements are present in sequences:
#' 
#' sqapply(sq_dna, function(sequence) sum(sequence == "A"))
#' sqapply(sq_ami, function(sequence) sum(sequence == "A"))
#' sqapply(sq_unt, function(sequence) sum(sequence == "A"))
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link[base]{lapply}}
#' @export
sqapply <- function(x, fun, ...,
                    single_string = FALSE, 
                    NA_letter = getOption("tidysq_NA_letter")) {
  assert_class(x, "sq")
  assert_function(fun)
  assert_flag(single_string)
  assert_string(NA_letter, min.chars = 1)
  
  CPP_apply_R_function(x, 
                       function(sequence) fun(sequence, ...), 
                       single_string, NA_letter)
}