#' Convert an object to sq
#' 
#' This generic function takes an object of arbitrary type and returns an \code{\link{sq}} object 
#' as an output. Default implementation of the method throws an error - there needs to be 
#' implemented a method for specified class in order for function to work.
#' 
#' @param x a object of a class that supports conversion to \code{\link{sq}}. 
#' @param ... other arguments passed to the method.
#' @return A \code{sq} object
#' 
#' @details In \pkg{tidysq} only for \code{\link{character}} vector there is an implemented method 
#' which in fact works exactly like \code{\link{construct_sq}} - you can also pass other arguments
#' like those supported by \code{\link{construct_sq}}.
#' 
#' @examples 
#' # constructing an example in usual way
#' sq_1 <- construct_sq("CTGA")
#' 
#' # using a method for character
#' sq_2 <- as.sq("CTGA")
#' 
#' # checking that both objects are identical
#' identical(sq_1, sq_2)
#' 
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
#' @export
as.sq <- function(x, ...) {
  UseMethod("as.sq")
}

#' @export
as.sq.default <- function(x, ...) {
  stop("'as.sq' cannot handle objects with this class")
}

#' @export
as.sq.character <- function(x, type = NULL, is_clean = NULL, non_standard = NULL, ...) {
  construct_sq(x, type, is_clean, non_standard)
}

#' Convert sq object into character vector
#' 
#' @description Coerce sequences from a \code{\link{sq}} object to \code{\link{character}} vector
#' of sequences
#' 
#' @param x a \code{\link{sq}} object to be converted
#' @param ... further arguments to be passed from or to other methods.
#' 
#' @return A \code{character} vector of the length the same as number
#' of sequences in the converted \code{sq} object
#' 
#' @details This method for class \code{\link{sq}} allows converting sequences from
#' the sq object into a character vector of length equal to the length 
#' of sq. Each element of resulting vector is a separate sequence. 
#' All attributes of the input sq are lost during the conversion to 
#' character vector.
#' 
#' @examples 
#' # Creating an object to work on:
#' sq_dna <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", 
#'                          "CAGACCANNNATAG"), type = 'dna')
#'                          
#' # Converting the sq object into a character vector:
#' as.character(sq_dna)
#'
#' @seealso sq
#' @export
as.character.sq <- function(x, ...) {
  vec_cast(x, character())
}

#' Convert sq object into matrix
#' 
#' @description Coerce sequences from a \code{\link{sq}} object to a 
#' \code{\link{matrix}}, in which rows correspond to sequences and columns to positions
#' 
#' @param x a \code{\link{sq}} object to be converted.
#' @param ... further arguments to be passed from or to other methods.
#' 
#' @return A \code{\link{matrix}} with number of rows the same as number of sequences
#' and number of columns corresponding to the length of the longest sequence
#' in the converted sq object. Matrix is either character (if type of \code{sq} is
#' \strong{ami}, \strong{dna}, \strong{rna}, \strong{atp} or \strong{unt})
#' or numeric (if type is \strong{enc}).
#' 
#' @details This method for class \code{sq} allows converting sequences from
#' the sq object into a matrix. Each row corresponds to the separate sequence
#' from the sq object, whereas each column indicates a single position within 
#' a sequence. Dimensions of matrix are determined by the number of sequences 
#' (rows) and the length of the longest sequence (columns). If a length of
#' sequence is smaller than the length of the longest sequence, the remaining
#' columns will be filled with \code{\link{NA}}. All attributes of the input
#' \code{sq} are lost during the conversion to matrix.
#' 
#' @examples 
#' # Creating objects to work on:
#' sq_same_len <- construct_sq(c("CGATAGACA", "TGACAAAAC", "GTGACCGTA"),
#'                             type = 'dna')
#' sq_diff_len <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", 
#'                               "CAGACCANNNATAG"), type = 'dna')
#' 
#' # Converting sq object containing sequences of the same lengths to matrix:
#' as.matrix(sq_same_len)
#' 
#' # Converting sq object containing sequences of different lengths to 
#' # matrix - NA introduced:
#' as.matrix(sq_diff_len)
#' 
#' @seealso \code{\link{sq}}
#' @export
as.matrix.sq <- function(x, ...) {
  max_len <- max(get_sq_lengths(x))
  ret <- do.call(rbind, lapply(.unpack_from_sq(x, "char"), function(row) row[1:max_len]))
  ret[ret == .get_na_char()] <- NA
  ret
}

#' @export
as.matrix.encsq <- function(x, ...) {
  ret <- NextMethod()
  storage.mode(ret) <- "numeric"
  ret
}

#' Check if object has specified type
#' 
#' Function to checks if object is a \code{\link{sq}} object without specifying type or
#' if it is a \code{\link{sq}} object with specific type.
#' @param x an object to be checked.
#' @return A \code{\link{logical}} value - \code{TRUE} if \code{x} has specified type, \code{FALSE}
#' otherwise.
#' 
#' @details 
#' These function does not only check objects classes - they also check if their format is 
#' correct (e.g. if they have alphabet parameter, if they have exactly one type, if they are
#' list of raws, etc. - to see details, how does \code{sq} object look like under the hood, 
#' read \code{\link[=sq]{sq class}} manual).
#' 
#' @examples 
#' sq_ami <- construct_sq(c("CVMPQGQQ", "AHLC--PPQ"))
#' sq_dna <- construct_sq(c("GGCAT", "TATC-A", "TGA"))
#' sq_rna <- construct_sq(c("CGAUUACG", "UUCUAGA", "UUCA"))
#' sq_unt <- construct_sq("BAHHAJJ&HAN&JD&")
#' sq_atp <- construct_sq(c("mALPVQAmAmA", "mAmAPQ"), non_standard = "mA")
#' sq_enc <- encode(sq_rna, c(A = 1.23, C = -0.72, G = 0.97, U = 3.01))
#' 
#' is.sq(sq_ami)
#' is.sq(sq_dna)
#' is.sq(sq_rna)
#' is.sq(sq_unt)
#' is.sq(sq_atp)
#' is.sq(sq_enc)
#' 
#' is.sq(c(1,2,3))
#' is.sq(LETTERS)
#' is.sq(TRUE)
#' is.sq(NULL)
#' 
#' is.amisq(sq_ami)
#' is.dnasq(sq_dna)
#' is.rnasq(sq_rna)
#' is.atpsq(sq_atp)
#' is.untsq(sq_unt)
#' is.encsq(sq_enc)
#' 
#' is.dnasq(sq_enc)
#' is.amisq(sq_atp)
#' is.untsq(sq_ami)
#' @seealso \code{\link{sq}}
#' @export
is.sq <- function(x) {
  tryCatch({.validate_sq(x); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @export
is.amisq <- function(x) {
  tryCatch({.validate_sq(x, type = "ami"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @export
is.dnasq <- function(x) {
  tryCatch({.validate_sq(x, type = "dna"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @export
is.rnasq <- function(x) {
  tryCatch({.validate_sq(x, type = "rna"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @export
is.untsq <- function(x) {
  tryCatch({.validate_sq(x, type = "unt"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @export
is.atpsq <- function(x) {
  tryCatch({.validate_sq(x, type = "atp"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @export
is.encsq <- function(x) {
  tryCatch({.validate_sq(x, type = "enc"); TRUE}, error = function(e) FALSE)
}

#' Compare sq objects
#' @description Compares input \code{\link{sq}} object with another given.
#' 
#' @details \code{`==`} converts left-hand side of the comparison (\code{x1}) to
#' a character vector using \code{\link{as.character}} and checks whether object
#' given on the right side can be compared with \code{\link{sq}} object. Function
#' also checks the type of \code{\link{sq}} object with which given object is to be
#' compared. If the type of \code{\link{sq}} object is "ami", "dna" or "rna" and
#' given sequence is a character vector with any lowercase letters, the function
#' replaces them with corresponding capital ones using \code{\link{toupper}}.
#' If right-hand side object (\code{x2}) is of class \code{\link{sq}}, then it is
#' converted to a character vector also using \code{\link{as.character}} function.
#' 
#' When both objects are already converted to character vectors, comparison is carried out 
#' element-wise with standard R rules, (e.g., recycling is used). You can check details at
#' \code{\link{Comparison}}.
#' 
#' Comparing sequences as characters vectors cause that various types of sequences
#' can be compared for example amino acid with nucleotide sequences so attention
#' should be paid, which sequence types are compared.
#' 
#' @param x1 an \code{\link{sq}} object.
#' @param x2 an object (a character vector or an \code{\link{sq}} object) to compare with
#' \code{x1}.
#' 
#' @return A \code{\link{logical}} vector indicating on which positions the objects are equal
#' 
#' @examples 
#' 
#' # Creating sq object to work on:
#' sq <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                      "CCCT", "CTGAATGT"), type = "dna")
#'                      
#' sq_different_len <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                                    "GGAA"), type = "dna")
#'                                    
#' sq_the_same_len <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                                  "CCCT", "CTGAATGT"), type = "dna")
#'                                                                                        
#' # Get an overview of the sequences:
#' summary(sq)
#' summary(sq_the_same_len)
#' summary(sq_different_len)
#'
#' # Comparing sq object with an object of the same length :
#' sq == sq_the_same_len
#' 
#' # Comparing object sq object with an object of a different length : 
#' sq == sq_different_len
#'  
#' # Comparing sq object to a given character vector of a different length:
#' sq == c('AAA','CCC')
#' 
#' # Comparing sq object to a given character vector of the same length:
#' sq == c("ACTGCTG", "CTTAGA",'CCCT', 'CTGAATGT')
#' 
#' # Comparing sq object to a given DNA element 'ATGTGA':
#' sq == 'ATGTGA'
#' 
#' # Comparing sq object to a given amino acid vector:
#' sq == c('RISGQQD','RISGQQD')
#'  
#' @seealso \code{\link{sq}} \code{\link{as.character}} \code{\link{is.sq}}
#' @export
`==.sq` <- function(x1, x2) {
  #TODO make it faster and lighter, maybe?
  if (is.sq(x2)) {
    x2 <- as.character(x2)
  } else if (!is.character(x2)) {
    stop("you cannot compare 'sq' object to object that is not character vector or 'sq' object")
  }
  
  if (.get_sq_type(x1) %in% c("ami", "dna", "rna")) {
    x2 <- toupper(x2)
  }
  
  as.character(x1) == x2
}

#' Get lengths of sequences in sq object
#' 
#' Function counts number of elements in each sequence in given \code{\link{sq}} object.
#' 
#' @param x an \code{\link{sq}} object.
#'  
#' @return A \code{\link{numeric}} vector, where each element gives length of according 
#' sequence from \code{\link{sq}} object.
#' 
#' @details This function allows returning numeric vector of lengths of sequences from
#' \code{\link{sq}} object. The numeric vector is as long as number of sequences present 
#' in \code{\link{sq}} object.
#' The function counts elements in all types of sequences.
#'
#' @examples 
#' # Creating an object to work on:
#' sq_dna <- construct_sq(c("ACGATTAGACG","GGATA"), type = "dna")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' 
#' # Counting number of elements in DNA sq object with defined type:
#' get_sq_lengths(sq_dna)
#' 
#' # Counting number of elements in amino acid sq object with defined type:
#' get_sq_lengths(sq_amino_acids)
#' 
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
#' @export
get_sq_lengths <- function(x) {
  if (length(x) == 0) numeric(0)
  else sapply(x, attr, "original_length")
}

#' Extract parts of a sq object
#' 
#' @description Operator to extract subsets of sq objects.
#' 
#' @param x sq object from which to extract element(s)
#' @param i,j,... indices specifying elements to extract. They may be 
#' \code{\link{numeric}}, \code{\link{character}} or \code{\link{logical}} vectors or empty. 
#' This function follows \code{\link[vctrs]{vctrs-package}} conventions regarding argument
#' interpretation for indexing vectors, which are a bit stricter that normal R
#' conventions, for example implicit argument recycling is prohibited.
#' 
#' @return \code{\link{sq}} object of the same type as input sq, containing
#' extracted elements
#' 
#' @details This function allows extracting specified sequences from the 
#' sq object and follows the normal R conventions. For details refer to the 
#' R documentation (see 
#' \url{https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Index-vectors}). 
#' Subsetting of the sq object does not affect its attributes (class and alphabet 
#' of the object). Attempt to extract elements using indices not present in
#' the object will return an error.
#' 
#' @examples 
#' # Creating objects to work on:
#' sq_unt <- construct_sq(c("AHSNLVSCTK$SH%&VS", "YQTVKA&#BSKJGY", 
#'                          "IAKVGDCTWCTY&GT", "AVYI#VSV&*DVGDJCFA"))
#' sq_ami <- construct_sq(c(s1 = "MAIATNCEPILLKNYAS", s2 = "YASDGLIPAKNTEWYTV", 
#'                          s3 = "TIKSNAALIETRY"), type = "ami")
#' 
#' # Subsetting using numeric vectors
#' # Extracting second element of the object:
#' sq_unt[2]
#' 
#' # Extracting elements from second to fourth:
#' sq_unt[2:4]
#' 
#' # Extracting all elements except the third:
#' sq_unt[-3]
#' 
#' # Extracting first and third element:
#' sq_unt[c(1,3)]
#' 
#' # Subsetting using character vectors
#' # Extracting elements named 's1' and 's3':
#' sq_ami[c('s1', 's3')]
#' 
#' # Subsetting using logical vectors
#' # Extracing first and third element:
#' sq_ami[c(TRUE, FALSE, TRUE)]
#' 
#' # Subsetting using empty vector
#' # Empty index will return all values:
#' sq_unt[]
#' 
#' @seealso \code{\link{sq}} \code{\link{bite}}
#' @name sqextract
#' @aliases sq-extract
NULL

#' Concatenate sq objects
#' 
#' @description Merges multiple \code{\link{sq}} and \code{character} objects into one larger
#' \code{sq}.
#' 
#' @param ... multiple \code{\link{sq}} and \code{character} objects. All of them have to have the same type and
#' subtype. If type is \strong{atp}, \strong{unt} or \strong{enc} also their alphabets have to
#' be exactly identical.
#' 
#' @return A \code{\link{sq}} object with length equal to sum of lengths of individual objects
#' passed as parameters. Elements of \code{\link{sq}} are concatenated just as if they were normal
#' lists (see \code{\link[base]{c}})
#' 
#' @examples 
#' sq_1 <- construct_sq(c("TAGACTAG", "", "CCGTAGATG"))
#' sq_2 <- construct_sq(c("TTGATAACG", "TGTATGTGA"))
#' sq_3 <- construct_sq(character(0))
#' sq_4 <- construct_sq("gaGG")
#' 
#' c(sq_1, sq_2, sq_3, sq_4)
#' 
#' @name sqconcatenate
#' @aliases sq-concatenate
NULL
