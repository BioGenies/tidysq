#' Subset sequences from sq objects
#' 
#' @description Extract a defined range of elements (amino acids or nucleotides) 
#' from a sequence.
#' 
#' @param sq \code{\link{sq}} object
#' @param indices \code{numeric} vector of subsequence indices to extract from
#' each sequence. The function follows the normal R conventions for indexing 
#' vectors, including negative indices.
#' 
#' @return \code{\link{sq}} object of the same type as input sq, where each 
#' element is a subsequence created by indexing corresponding sequence from 
#' input sq object with input indices.
#' 
#' @details Amino acids and nucleic acid sequences are represented as \code{\link{sq}} 
#' object in the \code{\link{tidysq}} package. Often one needs to get only a 
#' single letter or the sequence of a defined range from the original sequences. 
#' A subsequence is a sequence that can be derived from the original sequence 
#' by trimming some elements (letters) without changing the order of the 
#' remaining elements. To obtain a subsequence from each sequence contained in 
#' the \code{\link{sq}} object with the same indices. This is for example 
#' useful to extract an user-defined region from a sequence. 
#' 
#' The usage of \code{bite} follows the normal R conventions. For details 
#' refer to the R documentation (see 
#' \url{https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Index-vectors}). 
#' Extracting indices not present in the sequence results in introducing 
#' \code{\link{NA}} (‘Not Available’ / Missing Values) values. 
#' Information about it is printed on console depending on value of option 
#' 'tidysq_bite_na_action' - it can be either a warning (default), error, 
#' message or no information (you can check details in \code{\link{sq-options})}. 
#' \code{NA} values can be removed by using \code{\link{remove_na}} function.
#' 
#' @examples 
#' # Creating object, called sq to work on:
#' # The first four sequences of are miRNAs and isomiRs that were taken from 
#' # Tan et al. (2014) Nucleic Acid Research, doi: 10.1093/nar/gku656. 
#' # Selected miR-9 (isomiRs) sequences from neural progenitor stem cells
#' # look as follows: TCTTTGGTTATCTAGCTGTATGA, CTTTGGTTATCTAGCTGTATGA, 
#' # TCTTTGGTTATCTAGCTGTATG, TCTTTGGTTATCTAGCTGTATGAA
#' # The remaining sequences represent short artificial random sequences.
#' # Tan et al. have shown that "majority of miRNA genes encode mature isomers 
#' # that vary in size by one or more bases at the 3' and/or 5' end of the miRNA.
#' # In the following the bite function is used to trim nucleotieds 5' for a 
#' # simple visual alignment.
#' # construct_sq is use to create an sq object
#'
#' sq <- construct_sq(c("TCTTTGGTTATCTAGCTGTATGA", "CTTTGGTTATCTAGCTGTATGA", 
#'                      "TCTTTGGTTATCTAGCTGTATG", "TCTTTGGTTATCTAGCTGTATGAA", 
#'                      "ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), type = "nuc")
#'
#' # Get an overview of the sequences and show the first four only:
#' # The first four isomiRs sequences have lengths of 15, 14, 14 and 15
#' # nucleotides. The remaining artificial sequences have lengths of 5, 4, 3 
#' # and 5 nucleotides, respectively.
#' summary(sq)
#'
#' # Working with the first four miRNA/isomiR
#' # Removing first letter from the first four sequences (miRNAs/isomiRs).
#' # The sequence 1, 3 and 4 appear to be more similar.
#' bite(sq[1:4], -1)
#' 
#' # Extracting first five letters from each miRNA/ismomiR sequence:
#' bite(sq[1:4], 1:5)
#' 
#' # Working with all sequences
#' # extracting first letter from each sequence:
#' bite(sq, 1)
#' 
#' # Extracting first three letters from each sequence:
#' bite(sq, 1:3)
#' 
#' # Extracting second, fourth, third and second letters:
#' bite(sq, c(2,4,3,2))
#' 
#' # Extracting second to fifth letter - NA introduced:
#' bite(sq, 2:5)
#' 
#' # Extracting all from first to twentieth - NA introduced:
#' bite(sq, 1:20)
#' 
#' # Extracting all excluding first letter of sequence:
#' bite(sq, -1)
#' 
#' # Extracting all excluding second and sixth letter of sequence:
#' bite(sq, c(-2, -6))
#' 
#' 
#' @seealso sq remove_na sq-options
#' @export
bite <- function(sq, indices) {
  validate_sq(sq)

  .check_inds_are_numeric(indices)
  
  na_introduced <- FALSE
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  
  ret <- list(length(sq))
  for (i in 1:length(sq)) {
    s <- unpack_ints(sq[[i]], alph_size)
    s <- s[indices]
    if (any(is.na(s))) na_introduced <- TRUE
    s[is.na(s)] <- na_val
    ret[[i]] <- pack_ints(s, alph_size)
  }
  if (na_introduced) {
    .handle_opt_txt("tidysq_bite_na_action",
                    "some sequences are subsetted with index bigger than length - NA introduced")
  }
  .set_class_alph(ret, sq)
}
