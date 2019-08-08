#methods other than print of class sq

#' @exportMethod `[` sq
#' @export
`[.sq` <- function(x, i, j, ...) {
  ret <- NextMethod()
  ret <- lapply(ret, function(s) if (is.null(s)) raw(0) else s)
  class(ret) <- class(x)
  attr(ret, "alphabet") <- .get_alph(x)
  ret
}

#' @exportMethod as.character sq
#' @export
as.character.sq <- function(x, ...) {
  .debitify_sq(x, "string")
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

#' @exportMethod is sq
#' @export
is.sq <- function(x) {
  tryCatch({validate_sq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is amisq
#' @export
is.amisq <- function(x) {
  tryCatch({validate_amisq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is nucsq
#' @export
is.nucsq <- function(x) {
  tryCatch({validate_nucsq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is untsq
#' @export
is.untsq <- function(x) {
  tryCatch({validate_untsq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is atpsq
#' @export
is.atpsq <- function(x) {
  tryCatch({validate_atpsq(x); TRUE}, error = function(e) FALSE)
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
