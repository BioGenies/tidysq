#' Amino Acid Residue Background Frequencies
#' 
#' Amino acid frequency for the kingdoms of life in the Proteome-pI database with
#' added values for "All" from Proteome-pI database and from Seq2Logo method.
#' 
#' @name bg_freqs
#' @docType data
#' @format A data frame with 6 rows and 20 columns:
#' \describe{
#'   \item{name}{One-letter symbol for each of the 20 standard proteogenic amino acids}
#'   \item{Kingdom}{Viruses, Archaea, Bacteria, Eukaryota, All and Seq2logo}
#' }
#' @details 
#' Naturally observed amino acid residue background frequencies are available
#' from the Proteome Isoelectric Point Database. Proteome-pI Database is a database of
#' pre-computed isoelectric points for proteomes from different model organisms.
#' Values from Seq2Logo row are based on amino acid binding motifs and sequence profiles
#' including sequence weighting, pseudo counts and two-sided representation of 
#' amino acid enrichment and depletion.
#' 
#' 
#' @source Kozlowski LP. Proteome-pI: proteome isoelectric point database.
#' Nucleic Acids Res. 2017;45(D1):D1112–D1116. doi:10.1093/nar/gkw978
#' \url{https://academic.oup.com/nar/article/45/D1/D1112/2333931}
#' 
#' Thomsen MC, Nielsen M. Seq2Logo: a method for construction and visualization
#' of amino acid binding motifs and sequence profiles including sequence weighting,
#' pseudo counts and two-sided representation of amino acid enrichment and depletion.
#' Nucleic Acids Res. 2012;40(Web Server issue):W281–W287. doi:10.1093/nar/gks469
#' \url{https://academic.oup.com/nar/article/40/W1/W281/1076274}
#' 
#' @keywords datasets
#' @examples 
#' data(bg_freqs)
#' 
NULL
