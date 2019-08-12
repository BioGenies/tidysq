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
