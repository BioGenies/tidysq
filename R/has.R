#' Test sq object for presence of given motifs
#' 
#' @description Tests if elements of a \code{\link[=sq-class]{sq}} object
#' contain given motifs.
#' 
#' @template x
#' @param y [\code{character}]\cr
#'  Motifs to be searched for.
#' 
#' @return A \code{\link{logical}} vector of the same length as input \code{sq},
#' indicating which elements contain all given motifs.
#' 
#' @details
#' This function allows testing if elements of a \code{sq} object contain the
#' given motif or motifs. It returns a \code{logical} value for every element
#' of the \code{sq} object - \code{TRUE} if tested sequence contains searched
#' motif and \code{FALSE} otherwise. When multiple motifs are searched,
#' \code{TRUE} will be returned only for sequences that contain all given
#' motifs.
#' 
#' This function only indicates if a motif is present within a sequence, to find
#' all motifs and their positions within sequences use
#' \code{\link{find_motifs}}.
#' 
#' @template motif_details
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("ATGCAGGA", "GACCGNBAACGAN", "TGACGAGCTTAG"),
#'              alphabet = "dna_bsc")
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
#'              alphabet = c("mA", "mY", "nbA", "nsA"))
#'
#' # Testing if DNA sequences contain motif "ATG":
#' sq_dna %has% "ATG"
#'
#' # Testing if DNA sequences begin with "ATG":
#' sq_dna %has% "^ATG"
#'
#' # Testing if DNA sequences end with "TAG" (one of the stop codons):
#' sq_dna %has% "TAG$"
#'
#' # Test if amino acid sequences contain motif of two alanines followed by
#' # aspartic acid or asparagine ("AAB" motif matches "AAB", "AAD" and "AAN"):
#' sq_ami %has% "AAB"
#'
#' # Test if amino acid sequences contain both motifs:
#' sq_ami %has% c("AAXG", "MAT")
#'
#' # Test for sequences with multicharacter alphabet:
#' sq_atp %has% c("nsA", "mYmY$")
#'
#' @family bio_functions
#' @export
`%has%` <- function(x, y) {
  assert_character(y, any.missing = FALSE, min.len = 1)
  
  UseMethod("%has%")
}

#' @export
`%has%.default` <- function(x, y)
  stop_no_method(
    `%has%`, x,
    msg = function(cls) paste0("operator '%has%' is not overloaded for object of classes <",
                               paste0(class(x), collapse = ", "), ">")
  )

#' @export
`%has%.sq` <- function(x, y) {
  assert_alph_no_special_chars(alphabet(x))
  
  CPP_has(x, y, getOption("tidysq_NA_letter"))
}

#' @export
`%has%.sq_ami_ext` <- function(x, y) {
  y <- toupper(y)
  assert_motifs_for_type(y, "ami_ext")
  
  CPP_has(x, y, getOption("tidysq_NA_letter"))
}

#' @export
`%has%.sq_ami_bsc` <- `%has%.sq_ami_ext`

#' @export
`%has%.sq_dna_ext` <- function(x, y) {
  y <- toupper(y)
  assert_motifs_for_type(y, "dna_ext")
  
  CPP_has(x, y, getOption("tidysq_NA_letter"))
}

#' @export
`%has%.sq_dna_bsc` <- `%has%.sq_dna_ext`

#' @export
`%has%.sq_rna_ext` <- function(x, y) {
  y <- toupper(y)
  assert_motifs_for_type(y, "rna_ext")
  
  CPP_has(x, y, getOption("tidysq_NA_letter"))
}

#' @export
`%has%.sq_rna_bsc` <- `%has%.sq_rna_ext`
