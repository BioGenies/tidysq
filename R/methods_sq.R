#methods other than print of class sq
#' Extract parts of a sq object
#' 
#' @description Operator to extract subsets of sq objects.
#' 
#' @param x sq object from which to extract element(s)
#' @param i,j,... indices specifying elements to extract. They may be 
#' \code{numeric}, \code{character} or \code{logical} vectors or empty. 
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
#' seealso sq bite
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
as.sq.character <- function(x, type = NULL, is_clean = NULL, non_standard = NULL) {
  construct_sq(x, type, is_clean, non_standard)
}

#' @exportMethod as.character sq
#' @export
as.character.sq <- function(x, ...) {
  unlist(.debitify_sq(x, "string"))
}

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

#' @exportMethod is sq
#' @export
is.sq <- function(x) {
  tryCatch({validate_sq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is amisq
#' @export
is.amisq <- function(x) {
  tryCatch({validate_sq(x, type = "ami"); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is nucsq
#' @export
is.nucsq <- function(x) {
  tryCatch({validate_sq(x, type = "nuc"); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is untsq
#' @export
is.untsq <- function(x) {
  tryCatch({validate_sq(x, type = "unt"); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is atpsq
#' @export
is.atpsq <- function(x) {
  tryCatch({validate_sq(x, type = "atp"); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is encsq
#' @export
is.encsq <- function(x) {
  tryCatch({validate_sq(x, type = "enc"); TRUE}, error = function(e) FALSE)
}
#' Compare sq object 
#' @description Compares input \code{\link{sq}} object with another given.
#'   
#' @details \code{`==`} converts left hand side of comparision (x1) to chracters 
#' vector using \code{\link{as.character}} and checks whether given on the right side 
#' object can be compared with \code{\link{sq}} object. Function also check 
#' the type of \code{\link{sq}} object with which given object will be compared.
#' If the type of \code{\link{sq}} object is ami or nuc and given sequence  is 
#' character vector consisting lowercase, the function rewrites it into capital ones
#' with usage \code{\link{toupper}}. If right hand side object (x2) is \code{\link{sq}}
#' it is converted to character vector using also \code{\link{as.character}} function.
#' 
#' When both objects are already converted to character vectors, comparision is done 
#' elementwise with standard R rules, (e.g. recycling is used). You can check details 
#' \code{\link[Compare]{here}}.
#' 
#' Comparing sequences as characters vectors cause that various types of sequences
#' can be compared for example aminoacids with nucleotides sequences so attention 
#' should be paid which sequences types are compared. 
#' 
#' @param x1 \code{\link{sq}} object.
#' @param x2 an object (character vector or sq object) to compare with \code{\link{sq}}.
#' 
#' @return logical vector indicating on which positions objects are the same
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
#' @seealso sq as.character is.sq                                                          
#' @exportMethod `==` sq
#' @export
`==.sq` <- function(e1, e2) {
  #TODO make it faster and lighter, maybe?
  if (is.sq(e2)) {
    e2 <- as.character(e2)
  } else if (!is.character(e2)) {
    stop ("you cannot compare 'sq' object to object that is not character vector or 'sq' object")
  }
  
  type <- .get_sq_type(e1)
  if (type %in% c("ami", "nuc")) {
    e2 <- toupper(e2)
  }
  
  as.character(e1) == e2
}
