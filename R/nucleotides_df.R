#' Nucleotides abbreviations 
#' 
#' A dataset containing the nucleotide letter code and full names of nucleotides
#' according to IUPAC nomenclature.
#' @name nucleotides_df
#' @docType data 
#' @format  A data frame with 17 rows and 3 columns:
#' \describe{
#'  \item{one}{One-letter codes of nucleotides}
#'  \item{item}{Full name of the nucleotide}
#'  \item{amb}{Logical indicating if abbreviation is ambiguous, i.e., matches more than one nucleotide}
#'  }
#' @details The dataset contains a nucleotide alphabet of one-letter codes 
#' and full names of the nucleotides. It also includes a gap symbol '-' that can be found
#' in sequence alignments.
#' @source Johnson, A.D. (2010). An extended IUPAC nomenclature code for 
#' polymorphic nucleic acids. Bioinformatics 26, 1386–1389.
#' @keywords datasets
#' @examples 
#' data(nucleotides_df)
#'       
NULL          
