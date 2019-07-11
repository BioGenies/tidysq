#' Amino acid abbreviations
#' 
#' A dataset containing the one- and three-letter codes, and full names of amino acids according to the IUPAC nomenclature.
#' 
#' @name aminoacids_df
#' @docType data
#' @format A data frame with 27 rows and 4 columns:
#' \describe{
#'  \item{one}{One-letter codes of amino acids}
#'  \item{three}{Three-letter codes of amino acids}
#'  \item{full}{Full name of the amino acid}
#'  \item{amb}{Logical indicating if abbreviation is ambiguous, i.e., matches more than one amino acid}
#'  }
#' @source IUPAC
#' @keywords datasets
#' @examples 
#' data(aminoacids_df)
#' 
NULL