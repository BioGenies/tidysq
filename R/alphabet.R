#' Get alphabet of given sq object.
#' 
#' Function returns amino acid, DNA, RNA or atypical alphabet based on a \code{\link{sq}} 
#' object type. 
#' 
#' @param x a \code{\link{sq}} object to be recognized.
#'  
#' @return A character vector of letters of the alphabet.
#' 
#' @details This function allows returning alphabet of \code{sq} object which is a character
#' or numeric vector. The function reads provided \code{\link{sq}} object and determines,
#' which kind of sequences user assigned to a \code{\link{sq}} object (DNA, RNA, amino acid,
#' atypical or encoded one).
#' 
#' If \code{\link{sq}} has \strong{ami} type and \strong{cln} subtype, the function returns
#' the set of 20 aminoacids, the gap (-) and the stop codon (*) letter. If \strong{cln} is
#' missing, ambiguous letters are included as well.
#' 
#' If a \code{\link{sq}} contains \strong{dna} or \strong{rna} sequences and \strong{cln}
#' subtype, the function returns a set of respective 4 nucleotides with gap (-) element.
#' If \strong{cln} is missing, ambiguous letters are included as well.
#' 
#' If \code{\link{sq}} type is \strong{unt} or \strong{atp}, the function returns a list of 
#' letters present in sequences of a \code{\link{sq}} object.
#' 
#' If type is \strong{enc}, a numeric vector of values encoded for letters is returned
#' (see \code{\link{encode}}).
#' 
#' The details about amino acid and nucleotide alphabets can be checked in
#' \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}} respectively. General
#' information about alphabets and types of \code{sq} objects can be found in \code{\link{sq}}
#' class documentation.
#'
#' @seealso \code{\link{sq}} \code{\link{construct_sq}} \code{\link{encode}}
#' @export
alphabet <- function(x)
  attr(x, "alphabet")

`alphabet<-` <- function(x, value) {
  attr(x, "alphabet") <- value
  x
}

# alphabet creation ----

sq_alphabet <- function(alph, type) {
  # if exported add asserts
  new_vctr(
    alph,
    type = type,
    class = c("sq_alphabet", "character")
  )
}

sq_alphabet_ptype <- function(type)
  sq_alphabet(character(), type)

get_standard_alphabet <- function(type) {
  CPP_get_standard_alphabet(type)
}

obtain_alphabet <- function(x, sample_size = 4096, 
                            NA_letter = getOption("tidysq_NA_letter"),
                            ignore_case = FALSE) {
  CPP_obtain_alphabet(x, sample_size, NA_letter, ignore_case)
}

guess_standard_alphabet <- function(alph,
                                    NA_letter = getOption("tidysq_NA_letter")) {
  CPP_guess_standard_alph(alph, NA_letter)
}

# utility methods ----

`[.sq_alphabet` <- function(x, i,
                            NA_letter = getOption("tidysq_NA_letter")) {
  ret <- vec_data(x)[i]
  ret[i == (2 ^ size(x) - 1)] <- NA_letter
  ret
}

size <- function(alph) {
  ceiling(log2(length(alph) + 1))
}
