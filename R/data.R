#' Amino acid abbreviations
#' 
#' A dataset containing the one- and three-letter codes, and full names of 
#' amino acids according to the IUPAC nomenclature and HGVS Recommendations 
#' for the Description of Sequence Variants.
#' 
#' @name aminoacids_df
#' @docType data
#' @format A data frame with 27 rows and 4 columns:
#' \describe{
#'  \item{one}{One-letter codes of amino acids}
#'  \item{three}{Three-letter codes of amino acids}
#'  \item{full}{Full name of the amino acid}
#'  \item{amb}{Logical indicating if abbreviation is ambiguous, i.e., matches 
#'  more than one amino acid}
#'  }
#' @details
#' The dataset contains an amino acid alphabet of one-letter codes with the
#' corresponding three-letter codes and full names of the amino acids. It 
#' also includes a gap symbol '-' that can be found in sequence alignments
#' and termination symbol '*' used to indicate the end of translation.
#' 
#' @source (1984), Nomenclature and Symbolism for Amino Acids and Peptides. 
#' European Journal of Biochemistry 138, 9-37; Dunnen et al. (2016), HGVS 
#' Recommendations for the Description of Sequence Variants: 2016 Update. 
#' Human Mutation 37, 564-569.
#' 
#' @keywords datasets
#' @examples 
#' data(aminoacids_df)
#' 
NULL

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

#' The BLOSUM50 matrix
#' 
#' The BLOck SUbstition Matrix (BLOSUM) for sequences with less than 50\% similarity.
#' The matrix has been subset to the 20 proteogenic amino acids.
#' 
#' @name BLOSUM50
#' @docType data
#' @format A data frame with with 20 rows and 20 columns.
#' \describe{Contains an one-letter codes of amino acid as colums and rows names.
#' You can check three-letter codes and full names of the amino acids in
#' \code{\link{aminoacids_df}}.
#'  }
#' @details
#' BLOSUM matrix is a substitution matrix used for sequence alignment of proteins.
#' BLOSUM matrices are actual percentage identity values of sequences selected for
#' construction of the matrices. BLOSUM50 indicates that the sequences
#' selected for constructing the matrix share an average identity value of 50\%.
#' BLOSUM50 is good matrix for distantly related proteins.
#' The matrix made by matblas from blosum50.iij
#' BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
#' Cluster Percentage: >= 50
#' Entropy =   0.4808, Expected =  -0.3573
#' @source Henikoff, Steven, and Jorja G. Henikoff. 
#' "Amino acid substitution matrices from protein blocks." 
#' Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919.
#' \url{ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM50} 
#' 
#' @seealso aminoacids_df
#' @keywords datasets matrix
#' @examples 
#' data(BLOSUM50)
#' 
NULL

#' The BLOSUM62 matrix
#' 
#' The BLOck SUbstition Matrix (BLOSUM) for sequences with less than 62\% similarity.
#' The matrix has been subset to the 20 proteogenic amino acids.
#' 
#' @name BLOSUM62
#' @docType data
#' @format A data frame with with 20 rows and 20 columns.
#' \describe{Contains an one-letter codes of amino acid as colums and rows names.
#' You can check three-letter codes and full names of the amino acids in
#' \code{\link{aminoacids_df}}.
#'  }
#' @details
#' BLOSUM matrix is a substitution matrix used for sequence alignment of proteins.
#' BLOSUM matrices are actual percentage identity values of sequences selected for
#' construction of the matrices. BLOSUM62 indicates that the sequences
#' selected for constructing the matrix share an average identity value of 62\%.
#' BLOSUM62 is miderange matrix between close and  distangly related proteins.
#' Matrix made by matblas from blosum62.iij
#' BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#' Cluster Percentage: >= 62
#' Entropy = 0.6979, Expected = -0.5209
#' @source Henikoff, Steven, and Jorja G. Henikoff. 
#' "Amino acid substitution matrices from protein blocks." 
#' Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919.
#' \url{ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62} 
#' 
#' @seealso aminoacids_df
#' @keywords datasets matrix
#' @examples 
#' data(BLOSUM62)
#' 
NULL

#' The BLOSUM50_enc matrix
#' 
#' The BLOck SUbstition Matrix (BLOSUM) for sequences with less than 50\% similarity
#' transformed for encoding by dividing each value by 5.
#' The matrix has been subset to the 20 proteogenic amino acids.
#' 
#' @name BLOSUM50_enc
#' @docType data
#' @format A data frame with with 21 rows and 21 columns.
#' \describe{Contains an one-letter codes of amino acid as colums and rows names
#' and also row and column zero vectors representing 'X' (any amino acid).
#' You can check three-letter codes and full names of the amino acids in
#' \code{\link{aminoacids_df}}.
#'  }
#' @details
#' BLOSUM matrix is a substitution matrix used for sequence alignment of proteins.
#' BLOSUM matrices are actual percentage identity values of sequences selected for
#' construction of the matrices. BLOSUM50 indicates that the sequences
#' selected for constructing the matrix share an average identity value of 50\%.
#' BLOSUM50 is good matrix for distantly related proteins.
#' The matrix made by matblas from blosum50.iij
#' BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
#' Cluster Percentage: >= 50
#' Entropy =   0.4808, Expected =  -0.3573
#' @source Henikoff, Steven, and Jorja G. Henikoff. 
#' "Amino acid substitution matrices from protein blocks." 
#' Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919.
#' \url{ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM50} 
#' 
#' @seealso aminoacids_df Blosum50
#' @keywords datasets matrix
#' @examples 
#' data(BLOSUM50_enc)
#' 
NULL

#' The BLOSUM62_enc matrix
#' 
#' The BLOck SUbstition Matrix (BLOSUM) for sequences with less than 62\% similarity.
#' The matrix has been subset to the 20 proteogenic amino acids.
#' 
#' @name BLOSUM62_enc
#' @docType data
#' @format A data frame with with 21 rows and 21 columns.
#' \describe{Contains an one-letter codes of amino acid as colums and rows names
#' and also row and column zero vectors representing 'X' (any amino acid).
#' You can check three-letter codes and full names of the amino acids in
#' \code{\link{aminoacids_df}}.
#'  }
#' @details
#' BLOSUM matrix is a substitution matrix used for sequence alignment of proteins.
#' BLOSUM matrices are actual percentage identity values of sequences selected for
#' construction of the matrices. BLOSUM62 indicates that the sequences
#' selected for constructing the matrix share an average identity value of 62\%.
#' BLOSUM62 is miderange matrix between close and  distangly related proteins.
#' Matrix made by matblas from blosum62.iij
#' BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#' Cluster Percentage: >= 62
#' Entropy = 0.6979, Expected = -0.5209
#' @source Henikoff, Steven, and Jorja G. Henikoff. 
#' "Amino acid substitution matrices from protein blocks." 
#' Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919.
#' \url{ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62} 
#' 
#' @seealso aminoacids_df BLOSUM62
#' @keywords datasets matrix
#' @examples 
#' data(BLOSUM62_enc)
#' 
NULL

#' The BLOSUM50_pca matrix
#' 
#' 
#' This matrix enables principal component analysis (PCA) with usage \code{\link{BLOSUM50}}.
#' 
#' @name BLOSUM50_pca
#' @docType data
#' @format A data frame with with 21 rows and 20 columns.
#' \describe{Contains an one-letter codes of amino acid as rows (additionaly 'X' as any amino acid)
#' and ordinal numbers of principal component direction as columns names. Twenty columns represent 
#' the loadings of the twenty eigenvectors. 
#' You can check three-letter codes and full names of the amino acids in
#' \code{\link{aminoacids_df}}.
#'  }
#' @details
#' The BLOSUM50_pca matrix enables PCA calculation on proteins sequences aligments.
#' Components are generated by an eigenvector decomposition of the matrix formed
#' from pairwise similarity scores between each pair of sequences. The similarity score model
#' used for creating BLOSUM50_pca matrix is the \code{\link{BLOSUM50}}. 
#'  
#' @source Li, Jie & Koehl, Patrice. (2014). 3D representations of amino acids - 
#' Applications to protein sequence comparison and classification. 
#' Computational and Structural Biotechnology Journal. 11. 10.1016/j.csbj.2014.09.001. 
#' \url{ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM50} 
#' 
#' @seealso aminoacids_df Blosum50
#' @keywords datasets matrix
#' @examples 
#' data(BLOSUM50_pca)
#' 
NULL

#' The BLOSUM62_pca matrix
#' 
#' 
#' This matrix enables principal component analysis (PCA) with usage \code{\link{BLOSUM62}}.
#'
#' @name BLOSUM62_pca
#' @docType data
#' @format A data frame with with 21 rows and 20 columns.
#' \describe{Contains an one-letter codes of amino acid as rows (additionaly 'X' as any amino acid)
#' and ordinal numbers of principal component direction as columns names. Twenty columns represent 
#' the loadings of the twenty eigenvectors. 
#' You can check three-letter codes and full names of the amino acids in
#' \code{\link{aminoacids_df}}.
#'  }
#' @details
#' The BLOSUM62_pca matrix enables PCA calculation on proteins sequences aligments.
#' Components are generated by an eigenvector decomposition of the matrix formed
#' from pairwise similarity scores between each pair of sequences. The similarity score model
#' used for creating BLOSUM62_pca matrix is the \code{\link{BLOSUM62}}. 
#'  
#' @source Li, Jie & Koehl, Patrice. (2014). 3D representations of amino acids - 
#' Applications to protein sequence comparison and classification. 
#' Computational and Structural Biotechnology Journal. 11. 10.1016/j.csbj.2014.09.001. 
#' \url{ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62} 
#' 
#' @seealso aminoacids_df BLOSUM62
#' @keywords datasets matrix
#' @examples 
#' data(BLOSUM62_pca)
#' 
NULL

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

