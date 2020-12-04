#' Convert an object to sq
#' 
#' This generic function takes an object of arbitrary type and returns an \code{\link[=sq-class]{sq}} object
#' as an output. Default implementation of the method throws an error - there needs to be 
#' implemented a method for specified class in order for function to work.
#' 
#' @param x [\code{any}]
#'  An object of a class that supports conversion to \code{\link[=sq-class]{sq class}}.
#' @template three-dots
#' 
#' @return A \code{sq} object
#' 
#' @details In \pkg{tidysq} only for \code{\link{character}} vector there is an implemented method 
#' which in fact works exactly like \code{\link{construct_sq}} - you can also pass other arguments
#' like those supported by \code{\link{construct_sq}}.
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link{construct_sq}}
#' @export
# TODO: could we possibly delete this function?
as.sq <- function(x, ...)
  UseMethod("as.sq")

#' @export
as.sq.default <- function(x, ...)
  stop("'as.sq' cannot handle objects with this class")

#' @export
as.sq.character <- function(x, ...) sq(x, ...)

#' Convert sq object into character vector
#' 
#' @description Coerce sequences from a \code{\link[=sq-class]{sq}} object to \code{\link{character}} vector
#' of sequences
#' 
#' @template x
#' @template NA_letter
#' @template three-dots
#' 
#' @return A \code{character} vector of the length the same as number
#' of sequences in the converted \code{sq} object
#' 
#' @details This method for class \code{\link[=sq-class]{sq}} allows converting sequences from
#' the sq object into a character vector of length equal to the length 
#' of sq. Each element of resulting vector is a separate sequence. 
#' All attributes of the input sq are lost during the conversion to 
#' character vector.
#' 
#' @seealso sq
#' @export
as.character.sq <- function(x, ...,
                            NA_letter = getOption("tidysq_NA_letter"))
  unpack(x, "STRING", NA_letter)

#' Convert sq object into matrix
#' 
#' @description Coerce sequences from a \code{\link[=sq-class]{sq}} object to a
#' \code{\link{matrix}}, in which rows correspond to sequences and columns to positions
#' 
#' @template x
#' @template three-dots
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
#' @seealso \code{\link[=sq-class]{sq}}
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
#' Function to checks if object is a \code{\link[=sq-class]{sq}} object without specifying type or
#' if it is a \code{\link[=sq-class]{sq}} object with specific type.
#'
#' @template x
#'
#' @return A \code{\link{logical}} value - \code{TRUE} if \code{x} has specified type, \code{FALSE}
#' otherwise.
#' 
#' @details 
#' These function does not only check objects classes - they also check if their format is 
#' correct (e.g. if they have alphabet parameter, if they have exactly one type, if they are
#' list of raws, etc. - to see details, how does \code{sq} object look like under the hood, 
#' read \code{\link[=sq-class]{sq class}} manual).
#' 
#' @seealso \code{\link[=sq-class]{sq}}
#' @export
# TODO: could we possibly delete this function?
is.sq <- function(x)
  test_class(x, "sq")

# TODO: Are these necessary? Should we create e.g. is.sq_ami() that check for either one?
#' @rdname is.sq
#' @export
is.sq_ami_bsc <- function(x)
  test_class(x, "sq_ami_bsc")

#' Compare sq objects
#'
#' @description Compares input \code{\link[=sq-class]{sq}} object with either
#' another \code{sq} object or \code{character} vector.
#'
#' @param e1 [\code{sq}]\cr
#'  An object this comparison is applied to.
#' @param e2 [\code{sq} || \code{character}]\cr
#'  An object to compare with \code{x1}.
#' 
#' @return A \code{\link{logical}} vector indicating on which positions these
#' objects are equal.
#'
#' @details
#' \code{`==`} compares compatible object for equality of their respective
#' sequences. Objects are considered compatible, when either both have same
#' length or one of them is a scalar value (i.e. a vector of length 1).
#' Moreover, not every \code{e1} sq type can be compared to any \code{e2} sq
#' type.
#'
#' To see which types are compatible, see Details of
#' \code{\link{sq-concatenate}}.
#'
#' \code{`==`} returns logical vector, where each element describes whether
#' elements at position \code{n} of both \code{e1} and \code{e2} are equal in
#' meaning (that is, they may be represented differently, but their biological
#' interpretation must be indentical). If one of compared objects is a scalar,
#' then said logical vector describes comparison for each element of the other,
#' longer vector.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna_1 <- sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"),
#'                alphabet = "dna_bsc")
#' sq_dna_2 <- sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"),
#'                alphabet = "dna_bsc")
#' sq_dna_3 <- sq(c("ACTGCTG", "CTTAGA", "GGAA"),
#'                alphabet = "dna_bsc")
#' sq_dna_4 <- sq(c("ACTGCTG", "CTTAGA", "CCCT", "GTNANN"),
#'                alphabet = "dna_ext")
#' sq_ami_1 <- sq(c("ACTGCTG", "NIKAAR", "CCCT", "CTGAATGT"),
#'                alphabet = "ami_bsc")
#' sq_unt <- sq(c("AHSNLVSCTK$SH%&VS", "YQTVKA&#BSKJGY",
#'                "CCCT", "AVYI#VSV&*DVGDJCFA"))
#'
#' # Comparing sq object with an object of the same length:
#' sq_dna_1 == sq_dna_2
#' sq_dna_1 == c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT")
#'
#' # Cannot compare sq objects of different lengths:
#' \dontrun{
#' sq_dna_1 == sq_dna_3
#' sq_dna_1 == c("AAA", "CCC")
#' }
#'
#' # Unless comparing sq object with scalar value:
#' sq_dna_1 == "CTTAGA"
#'
#' # It's possible to compare basic and extended types:
#' sq_dna_1 == sq_dna_4
#'
#' # Mixing DNA, RNA and amino acid types throws an error, however:
#' \dontrun{
#' sq_dna_1 == sq_ami_1
#' }
#'
#' # On the other hand, unt sq is acceptable everywhere:
#' sq_dna_1 == sq_unt
#' sq_dna_4 == sq_unt
#' sq_ami_1 == sq_unt
#'
#' @family util_functions
#' @export
`==.sq` <- function(e1, e2) vec_equal(e1, e2)

#' Get lengths of sequences in sq object
#' 
#' @description Returns number of elements in each sequence in given
#' \code{\link[=sq-class]{sq}} object.
#' 
#' @template x
#'  
#' @return A \code{\link{numeric}} vector, where each element gives length of
#' corresponding sequence from \code{\link[=sq-class]{sq}} object.
#' 
#' @details
#' Due to storage implementation, using \code{\link[base]{lengths}} method
#' returns length of stored raw vectors instead of real sequence lengths. This
#' function accesses \code{original_length} attribute of each sequence, which
#' attribute stores information about how many elements are there in given
#' sequence.
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
#'              alphabet = "dna_ext")
#'
#' # Counting number of elements in sq object:
#' get_sq_lengths(sq_dna)
#' get_sq_lengths(sq_ami)
#'
#' @family util_functions
#' @export
get_sq_lengths <- function(x) {
  if (length(x) == 0) numeric(0)
  else vapply(x, attr, numeric(1), "original_length")
}

#' Extract parts of a sq object
#' 
#' @description Operator to extract subsets of sq objects.
#' 
#' @template x
#' @param i,j,... [\code{numeric} || \code{logical}]\cr
#'  Indices specifying elements to extract.
#' 
#' @return \code{\link[=sq-class]{sq}} object of the same type as the input,
#' containing extracted elements
#' 
#' @details
#' This function follows \code{\link[vctrs]{vctrs-package}} conventions
#' regarding argument interpretation for indexing vectors, which are a bit
#' stricter that normal R conventions, for example implicit argument recycling
#' is prohibited. Subsetting of the \code{sq} object does not affect its
#' attributes (class and alphabet of the object). Attempt to extract elements
#' using indices not present in the object will return an error.
#'
#' @examples
#' # Creating object to work on:
#' sq_unt <- sq(c("AHSNLVSCTK$SH%&VS", "YQTVKA&#BSKJGY",
#'                "IAKVGDCTWCTY&GT", "AVYI#VSV&*DVGDJCFA"))
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
#' # Subsetting using logical vectors
#' # Extracing first and third element:
#' sq_unt[c(TRUE, FALSE, TRUE, FALSE)]
#'
#' # Subsetting using empty vector returns all values:
#' sq_unt[]
#'
#' # Using NULL, on the other hand, returns empty sq:
#' sq_unt[NULL]
#'
#' @family util_functions
#' @name sqextract
#' @aliases sq-extract
NULL

#' Concatenate sq objects
#' 
#' @description Merges multiple \code{\link[=sq-class]{sq}} and possibly
#' \code{character} objects into one larger \code{sq} object.
#' 
#' @param ... [\code{sq} || \code{character}]\cr
#'  Multiple objects. For exact behaviour, check Details section. First argument
#'  must be of \code{sq} class due to R mechanism of single dispatch. If this is
#'  a problem, recommended alternative is \code{\link[vctrs]{vec_c}} method from
#'  \code{\link[vctrs]{vctrs-package}} package.
#' 
#' @return \code{\link[=sq-class]{sq}} object with length equal to sum of
#' lengths of individual objects passed as parameters. Elements of
#' \code{\link[=sq-class]{sq}} are concatenated just as if they were normal
#' lists (see \code{\link[base]{c}}).
#' 
#' @details
#' Whenever all passed objects are of one of standard types (that is,
#' \strong{dna_bsc}, \strong{dna_ext}, \strong{rna_bsc}, \strong{rna_ext},
#' \strong{ami_bsc} or \strong{ami_ext}), returned object is of the same class,
#' as no changes to alphabet are needed.
#'
#' It's possible to mix both basic and extended types within one call to
#' \code{c()}, however they all must be of the same type (that is, either
#' \strong{dna}, \strong{rna} or \strong{ami}). In this case, returned object
#' is of extended type.
#'
#' Mixing \strong{dna}, \strong{rna} and \strong{ami} types is prohibited, as
#' interpretation of letters differ depending on the type.
#'
#' Whenever all objects are either of \strong{atp} type, returned object is also
#' of this class and resulting alphabet is equal to set union of all input
#' alphabets.
#'
#' \strong{unt} type can be mixed with any other type, resulting in \strong{unt}
#' object with alphabet equal to set union of all input alphabets. In this case,
#' it is possible to concatenate \strong{dna} and \strong{ami} objects, for
#' instance, by concatenating one of them first with \strong{unt} object.
#' However, it is strongly discouraged, as it may result in unwanted
#' concatenation of DNA and amino acid sequences.
#'
#' Whenever a character vector appears, it does not influence resulting sq type.
#' Each element is treated as separate sequence. If any of letters in this
#' vector does not appear in resulting alphabet, it is silently replaced with
#' \code{NA}.
#'
#' Due to R dispatch mechanism passing character vector as first will return
#' class-less list. This behaviour is effectively impossible and definitely
#' unrecommended to fix, as fixing it would involve changing \code{c} primitive.
#' If such possibility is necessary, \code{\link[vctrs]{vec_c}} is a better
#' alternative.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna_1 <- sq(c("GGACTGCA", "CTAGTA", ""), alphabet = "dna_bsc")
#' sq_dna_2 <- sq(c("ATGACA", "AC-G", "-CCAT"), alphabet = "dna_bsc")
#' sq_dna_3 <- sq(character(), alphabet = "dna_bsc")
#' sq_dna_4 <- sq(c("BNACV", "GDBADHH"), alphabet = "dna_ext")
#' sq_rna_1 <- sq(c("UAUGCA", "UAGCCG"), alphabet = "rna_bsc")
#' sq_rna_2 <- sq(c("-AHVRYA", "G-U-HYR"), alphabet = "rna_ext")
#' sq_rna_3 <- sq("AUHUCHYRBNN--", alphabet = "rna_ext")
#' sq_ami <- sq("ACHNK-IFK-VYW", alphabet = "ami_bsc")
#' sq_unt <- sq("AF:gf;PPQ^&XN")
#'
#' # Concatenating dna_bsc sequences:
#' c(sq_dna_1, sq_dna_2, sq_dna_3)
#' # Concatenating rna_ext sequences:
#' c(sq_rna_2, sq_rna_3)
#' # Mixing dna_bsc and dna_ext:
#' c(sq_dna_1, sq_dna_4, sq_dna_2)
#'
#' # Mixing DNA and RNA sequences doesn't work:
#' \dontrun{
#' c(sq_dna_3, sq_rna_1)
#' }
#'
#' # untsq can be mixed with DNA, RNA and amino acids:
#' c(sq_ami, sq_unt)
#' c(sq_unt, sq_rna_1, sq_rna_2)
#' c(sq_dna_2, sq_unt, sq_dna_3)
#'
#' # Character vectors are also acceptable:
#' c(sq_dna_2, "TGCA-GA")
#' c(sq_rna_2, c("UACUGGGACUG", "AUGUBNAABNRYYRAU"), sq_rna_3)
#' c(sq_unt, "&#JIA$O02t30,9ec", sq_ami)
#'
#' @family util_functions
#' @name sqconcatenate
#' @aliases sq-concatenate
NULL
