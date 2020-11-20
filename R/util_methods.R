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
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
#' @export
as.sq <- function(x, ...)
  UseMethod("as.sq")

#' @export
as.sq.default <- function(x, ...)
  stop("'as.sq' cannot handle objects with this class")

#' @export
as.sq.character <- function(x,
                            alphabet = guess_sq_type(x),
                            NA_letter = getOption("tidysq_NA_letter"),
                            safe_mode = getOption("tidysq_safe_mode"))
  sq(x, alphabet, NA_letter, safe_mode)

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
#' @seealso sq
#' @export
as.character.sq <- function(x, ...)
  vec_cast(x, character())

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
#' @seealso \code{\link{sq}}
#' @export
as.matrix.sq <- function(x, ...) {
  max_len <- max(get_sq_lengths(x))
  ret <- do.call(rbind, lapply(unpack(x, "STRINGS"), function(row) row[1:max_len]))
  ret[ret == getOption("tidysq_NA_letter")] <- NA
  ret
}

#' @export
as.matrix.sq_enc <- function(x, ...) {
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
#' @seealso \code{\link{sq}}
#' @export
is.sq <- function(x)
  test_class(x, "sq")

# TODO: Are these necessary? Should we create e.g. is.sq_ami() that check for either one?
#' @rdname is.sq
#' @export
is.sq_ami_bsc <- function(x)
  test_class(x, "sq_ami_bsc")

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
#' @seealso \code{\link{sq}} \code{\link{as.character}} \code{\link{is.sq}}
#' @export
`==.sq` <- vec_equal

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
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
#' @export
get_sq_lengths <- function(x) {
  if (length(x) == 0) numeric(0)
  else vapply(x, attr, numeric(1), "original_length")
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
#' @seealso \code{\link{sq}} \code{\link{bite}}
#' @name sqextract
#' @aliases sq-extract
NULL

#' Concatenate sq objects
#' 
#' @description Merges multiple \code{\link{sq}} and maybe \code{character} objects
#' into one larger \code{sq} object.
#' 
#' @param ... multiple \code{\link{sq}} and \code{character} objects. For exact behaviour,
#' check Details section. First argument must be of \code{sq} class due to R mechanism of
#' single dispatch. If this is a problem, recommended alternative is \code{\link[vctrs]{vec_c}}
#' method from \code{vctrs} package.
#' 
#' @return A \code{\link{sq}} object with length equal to sum of lengths of individual objects
#' passed as parameters. Elements of \code{\link{sq}} are concatenated just as if they were normal
#' lists (see \code{\link[base]{c}})
#' 
#' @details Whenever all passed objects are either \code{dnasq}, \code{rnasq} or \code{amisq},
#' returned object is also of the same class. If all of them contain \code{clnsq} subtype,
#' resulting object contains it as well, else \code{clnsq} subtype is dropped.
#' 
#' Mixing \code{dnasq}, \code{rnasq} and \code{amisq} is prohibited, as interpretation of
#' symbols differ depending on the type.
#' 
#' Whenever all objects are either \code{untsq} or \code{atpsq}, returned object is also
#' of the same class. These classes are not used with \code{clnsq} subtype, so it never
#' appears in this context. Resulting object has alphabet equal to the set union of all
#' alphabets of the involved objects.
#' 
#' Moreover, \code{untsq} objects may be concatenated with \code{dnasq}, \code{rnasq}
#' and \code{amisq} objects, resulting in an \code{untsq} object with alphabet equal
#' to the set union of all alphabets involved. However, user is strongly encouraged
#' to use this possibility with caution, as it may result in unwanted concatenation
#' of DNA and amino acid sequences.
#' 
#' Whenever character vectors are passed as second and later argument, they don't influence
#' resulting \code{sq} object type. Each element of vector is interpreted as separate
#' sequence. If resulting \code{sq} is predicted to have \code{clnsq} subtype, then passing
#' character vector with characters not included in the resulting alphabet will cause code
#' to fail with an exception. If result doesn't have \code{clnsq} subtype, all "foreign"
#' characters will be silently replaced will \code{NA}.
#' 
#' Due to R dispatch mechanism passing character vector as first will return class-less
#' list. This behaviour is effectively impossible and definitely unrecommended to fix, as
#' fixing it would involve changing \code{c} primitive. If such possibility is necessary,
#' \code{\link[vctrs]{vec_c}} is a better alternative.
#'
#' @name sqconcatenate
#' @aliases sq-concatenate
NULL
