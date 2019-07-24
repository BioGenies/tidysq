#' Subset sequences from sq
#' 
#' Obtain subsequence from each sequence contained in the sq object with the 
#' same indices.
#' 
#' @param sq \code{\link{sq}} object
#' @param indices \code{numeric} vector of subsequence indices to extract from
#' each sequence. Follows the normal R conventions for indexing vectors, 
#' including negative indices.
#' 
#' @return \code{\link{sq}} object of the same type as input sq, where each 
#' element is a subsequence created by indexing corresponding sequence from 
#' input sq object with input indices.
#' 
#' @details This function follows the normal R conventions, thus extracting 
#' indices not present in the sequence results in introducing NA values. 
#' Information about it is printed on console depending on value of option 
#' 'tidysq_bite_na_action' - it can be either a warning (default), error, 
#' message or no information (you can check details in \link{sq-options}). 
#' NA values can be removed by using \link{remove_na} function.
#' 
#' @examples 
#' # creating object to work on:
#' sq <- construct_sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), type = "nuc")
#' 
#' # extracting first letter from each sequence:
#' bite(sq, 1)
#' 
#' # extracting first three letters from each sequence:
#' bite(sq, 1:3)
#' 
#' # extracting second, fourth, third and second letters:
#' bite(sq, c(2,4,3,2))
#' 
#' # extracting second to fifth letter - NA introduced:
#' bite(sq, 2:5)
#' 
#' # extracting all from first to twentieth - NA introduced:
#' bite(sq, 1:20)
#' 
#' # extracting all excluding first letter of sequence:
#' bite(sq, -1)
#' 
#' # extracting all excluding second and sixth letter of sequence:
#' bite(sq, c(-2, -6))
#' 
#' 
#' @seealso sq remove_na sq-options
#' @export
bite <- function(sq, indices) {
  validate_sq(sq)

  if (!(is.numeric(indices) && 
        floor(indices) == indices)) {
    stop("'indices' has to be an integer vector")
  }
  
  na_introduced <- FALSE
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  
  ret <- list(length(sq))
  for (i in 1:length(sq)) {
    s <- .bit_to_int(sq[[i]], alph_size)
    s <- s[indices]
    if (any(is.na(s))) na_introduced <- TRUE
    s[is.na(s)] <- na_val
    ret[[i]] <- .int_to_bit(s, alph_size)
  }
  if (na_introduced) {
    .handle_opt_txt("tidysq_bite_na_action",
                    "some sequences are subsetted with index bigger than length - NA introduced")
  }
  .set_class_alph(ret, sq)
}
