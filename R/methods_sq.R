#' Extract parts of a sq object
#' 
#' @rdname sqextract
#' @aliases sq-extract
#' @description Operator to extract subsets of sq objects.
#' 
#' @param x sq object from which to extract element(s)
#' @param i,j,... indices specifying elements to extract. They may be 
#' \code{\link{numeric}}, \code{\link{character}} or \code{\link{logical}} vectors or empty. 
#' This function follows normal R conventions for indexing vectors, including 
#' negative indices.
#' 
#' @return \code{\link{sq}} object of the same type as input sq, containing
#' extracted elements
#' 
#' @details This function allows extracting specified sequences from the 
#' sq object and follows the normal R conventions. For details refer to the 
#' R documentation (see 
#' \url{https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Index-vectors}). 
#' Subsetting of the sq object will not drop its attributes (class and alphabet 
#' of the object). Numeric values are coerced to integer as by 
#' \code{\link{as.integer}} and hence truncated towards zero. In case of 
#' logical vectors indicating elements/slices to select, they are recycled if 
#' necessary to match the length of the sq object. Attempt to extract elements 
#' using indices not present in the object will return an error.  
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
#' # Extracting using non-integer indices - truncation:
#' i <- 2.754
#' sq_unt[i]
#' 
#' # Subsetting using character vectors
#' # Extracting elements named 's1' and 's3':
#' sq_ami[c('s1', 's3')]
#' 
#' # Subsetting using logical vectors
#' # Extracing first and third element:
#' sq_ami[c(TRUE, FALSE, TRUE)]
#'
#' # Extracting every other element - vector will be recycled:
#' sq_unt[c(FALSE, TRUE)]
#' 
#' # Subsetting using empty vector
#' # Empty index will return all values:
#' sq_unt[]
#' 
#' @seealso \code{\link{sq}} \code{\link{bite}}
#'   
#' @exportMethod `[` sq
#' @export
`[.sq` <- function(x, i, j, ...) {
  ret <- NextMethod()
  ret <- lapply(ret, function(s) if (is.null(s)) raw(0) else s)
  class(ret) <- class(x)
  attr(ret, "alphabet") <- .get_alph(x)
  ret
}

#' Concatenate sq objects
#' 
#' Merges multiple \code{\link{sq}} objects of the same type into one larger.
#' 
#' @param ... multiple \code{\link{sq}} objects. All of them have to have the same type and
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
#' @exportMethod c sq
#' @export
c.sq <- function(...) {
  args <- list(...)
  if (!all(sapply(args, is.sq)))
    stop("not all elements passed to function are 'sq' objects")
  types <- sapply(args, .get_sq_type)
  if (length(unique(unlist(types))) != 1)
    stop("not all sq objects have the same type")
  ret <- unlist(args, recursive = FALSE)
  alphs <- lapply(args, .get_alph)
  if (!all(sapply(alphs[-1], function(alph) identical(alphs[[1]], alph))))
    stop("all of sq objects should have identical alphabets")
  .set_class_alph(ret, args[[1]])
}

#' Convert an object to sq
#' 
#' This generic function takes an object of arbitrary type and returns a \code{\link{sq}} object 
#' as an output. Default implementation of the method throws an error - there needs to be 
#' implemented a method for specified class in order for function to work.
#' 
#' @param x a object of a class that supports conversion to \code{\link{sq}}. 
#' @param ... other arguments passed to the method.
#' @return A \code{sq} object
#' 
#' @details In \pkg{tidysq} only for \code{\link{character}} vector there is implemented a method 
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

#' @exportMethod as.sq default
#' @export
as.sq.default <- function(x, ...) {
  stop("'as.sq' cannot handle objects with this class")
}

#' @exportMethod as.sq character
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
#' sq_nuc <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", 
#'                          "CAGACCANNNATAG"), type = 'nuc')
#'                          
#' # Converting the sq object into a character vector:
#' as.character(sq_nuc)
#'
#' @seealso sq
#' @exportMethod as.character sq
#' @export
as.character.sq <- function(x, ...) {
  unlist(.debitify_sq(x, "string"))
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
#' \strong{ami}, \strong{nuc}, \strong{atp} or \strong{unt}) or numeric (if type
#' is \strong{enc}).
#' 
#' @details This method for class \code{sq} allows converting sequences from
#' the sq object into a matrix. Each row corresponds to the separate sequence
#' from the sq object, whereas each column indicates a single position within 
#' a sequence. Dimensions of matrix are determined by the number of sequences 
#' (rows) and the length of the longest sequence (columns). If a length of
#' sequence is smaller than the length of the longest sequence, the remaining
#' columns will be filled with \code{\link{NA}}. All attributes of the input \code{sq} are lost 
#' during the conversion to matrix.
#' 
#' @examples 
#' # Creating objects to work on:
#' sq_same_len <- construct_sq(c("CGATAGACA", "TGACAAAAC", "GTGACCGTA"),
#'                             type = 'nuc')
#' sq_diff_len <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", 
#'                               "CAGACCANNNATAG"), type = 'nuc')
#' 
#' # Converting sq object containing sequences of the same lengths to matrix:
#' as.matrix(sq_same_len)
#' 
#' # Converting sq object containing sequences of different lengths to 
#' # matrix - NA introduced:
#' as.matrix(sq_diff_len)
#' 
#' @seealso \code{\link{sq}}
#' @exportMethod as.matrix sq
#' @export
as.matrix.sq <- function(x, ...) {
  x <- .debitify_sq(x, "char")
  max_len <- max(lengths(x))
  ret <- do.call(rbind, lapply(x, function(row) row[1:max_len]))
  ret[ret == .get_na_char()] <- NA
  ret
}


#' @exportMethod as.matrix encsq
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
#' sq_nuc <- construct_sq(c("CGAUUACG", "UUCUAGA", "UUCA"))
#' sq_unt <- construct_sq("BAHHAJJ&HAN&JD&")
#' sq_atp <- construct_sq(c("mALPVQAmAmA", "mAmAPQ"), non_standard = "mA")
#' sq_enc <- encode(sq_nuc, c(A = 1.23, C = -0.72, G = 0.97, U = 3.01))
#' 
#' is.sq(sq_ami)
#' is.sq(sq_nuc)
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
#' is.nucsq(sq_nuc)
#' is.atpsq(sq_atp)
#' is.untsq(sq_unt)
#' is.encsq(sq_enc)
#' 
#' is.nucsq(sq_enc)
#' is.amisq(sq_atp)
#' is.untsq(sq_ami)
#' @seealso \code{\link{sq}}
#' @exportMethod is sq
#' @export
is.sq <- function(x) {
  tryCatch({validate_sq(x); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @exportMethod is amisq
#' @export
is.amisq <- function(x) {
  tryCatch({validate_sq(x, type = "ami"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @exportMethod is nucsq
#' @export
is.nucsq <- function(x) {
  tryCatch({validate_sq(x, type = "nuc"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @exportMethod is untsq
#' @export
is.untsq <- function(x) {
  tryCatch({validate_sq(x, type = "unt"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @exportMethod is atpsq
#' @export
is.atpsq <- function(x) {
  tryCatch({validate_sq(x, type = "atp"); TRUE}, error = function(e) FALSE)
}

#' @rdname is.sq
#' @exportMethod is encsq
#' @export
is.encsq <- function(x) {
  tryCatch({validate_sq(x, type = "enc"); TRUE}, error = function(e) FALSE)
}
#' Compare sq object 
#' @description Compares input \code{\link{sq}} object with another given.
#'   
#' @details \code{`==`} converts left hand side of comparison (\code{x1}) to characters 
#' vector using \code{\link{as.character}} and checks whether given on the right side 
#' object can be compared with \code{\link{sq}} object. Function also check 
#' the type of \code{\link{sq}} object with which given object will be compared.
#' If the type of \code{\link{sq}} object is ami or nuc and given sequence is 
#' character vector consisting lowercase, the function rewrites it into capital ones
#' with usage \code{\link{toupper}}. If right-hand side object (x2) is \code{\link{sq}}
#' it is converted to character vector using also \code{\link{as.character}} function.
#' 
#' When both objects are already converted to character vectors, comparison is carried out 
#' element-wise with standard R rules, (e.g., recycling is used). You can check details 
#' \code{\link[Compare]{here}}.
#' 
#' Comparing sequences as characters vectors cause that various types of sequences
#' can be compared for example amino acids with nucleotides sequences so attention 
#' should be paid, which sequences types are compared. 
#' 
#' @param x1 a \code{\link{sq}} object.
#' @param x2 an object (character vector or sq object) to compare with \code{\link{sq}}.
#' 
#' @return A \code{\link{logical}} vector indicating on which positions objects are the same
#' 
#' @examples 
#' 
#' # Creating sq object to work on:
#' sq <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                      "CCCT", "CTGAATGT"), type = "nuc")
#'                      
#' sq_different_len <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                                    "GGAA"), type = "nuc")
#'                                    
#' sq_the_same_len <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                                  "CCCT", "CTGAATGT"), type = "nuc")
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
#' # Comparing sq object to given character vector of a different length:
#' sq == c('AAA','CCC')
#' 
#' # Comparing sq object to given character vector of a the same length:
#' sq == c("ACTGCTG", "CTTAGA",'CCCT', 'CTGAATGT')
#' 
#' # Comparing sq object to given nucleotide element 'ATGTGA':
#' sq == 'ATGTGA'
#' 
#' # Comparing sq object to given amino acids vector:
#' sq == c('RISGQQD','RISGQQD')
#'  
#' @seealso \code{\link{sq}} \code{\link{as.character}} \code{\link{is.sq}}                                                         
#' @exportMethod `==` sq
#' @export
`==.sq` <- function(x1, x2) {
  #TODO make it faster and lighter, maybe?
  if (is.sq(x2)) {
    x2 <- as.character(x2)
  } else if (!is.character(x2)) {
    stop("you cannot compare 'sq' object to object that is not character vector or 'sq' object")
  }
  
  type <- .get_sq_type(x1)
  if (type %in% c("ami", "nuc")) {
    x2 <- toupper(x2)
  }
  
  as.character(x1) == x2
}
