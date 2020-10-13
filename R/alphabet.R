# alphabet assignment ----

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
#' @examples 
#' # Creating an object to work on:
#' sq_dna <- construct_sq(c("ACGATTAGACG","GGATA"), type = "dna")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_untyped <- construct_sq(c("ACGA&&TTAGACG&"), type = "unt")
#' 
#' # Testing DNA sq object with defined type:
#' get_sq_alphabet(sq_dna)
#' 
#' # Testing amino acid sq object with defined type:
#' get_sq_alphabet(sq_amino_acids)
#' 
#' # Testing nucleotide sq object without defined type:
#' get_sq_alphabet(sq_untyped)
#' 
#'   
#' @seealso \code{\link{sq}} \code{\link{construct_sq}} \code{\link{encode}}
#' @export
get_sq_alphabet <- function(x) {
  assert_class(x, "sq")
  alphabet(x)
}

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

sq_alphabet_ptype <- function()
  sq_alphabet(character(), character())

get_standard_alph <- function(type) {
  sq_alphabet(
    switch (type,
            dna_bsc = nucleotides_df[nucleotides_df[["dna"]], "one"],
            dna_ext = nucleotides_df[nucleotides_df[["dna"]] | nucleotides_df[["amb"]], "one"],
            rna_bsc = nucleotides_df[nucleotides_df[["rna"]], "one"],
            rna_ext = nucleotides_df[nucleotides_df[["rna"]] | nucleotides_df[["amb"]], "one"],
            ami_bsc = aminoacids_df[!aminoacids_df[["amb"]], "one"],
            ami_ext = aminoacids_df[, "one"]
    ),
    type
  )
}

.skip_characters <- function(alph, chars)
  vec_restore(setdiff(alph, chars), alph)

# alphabet reading ----

`[.sq_alphabet` <- function(x, i,
                            NA_letter = getOption("tidysq_NA_letter")) {
  ret <- vec_data(x)[i]
  ret[i == .get_na_val(x)] <- NA_letter
  ret
}

# various internal methods put together (to check!) ----

.get_alph_size <- function(alph) {
  ceiling(log2(length(alph) + 2))
}

.get_na_val <- function(alph) {
  2 ^ .get_alph_size(alph) - 1
}

.get_real_alph <- function(str_sq) {
  new_vctr(
    C_get_real_alph(str_sq),
    na_letter = getOption("tidysq_NA_letter"),
    class = c("sq_alphabet", "character")
  )
}
