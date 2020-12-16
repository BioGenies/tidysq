#' Convert an object to sq
#' 
#' @description Takes an object of arbitrary type and returns an
#' \code{\link[=sq-class]{sq}} object as an output.
#' 
#' @param x [\code{any}]
#'  An object of a class that supports conversion to \code{sq} class.
#' @template three-dots
#' 
#' @return An \code{sq} object.
#' 
#' @details
#' There are two possible cases: if \code{x} is a character vector, then this
#' method calls \code{\link{sq}} function, else it passes \code{x} to
#' \code{\link{import_sq}} and hopes it works.
#'
#' @examples
#' # Constructing an example sequence in the usual way:
#' sq_1 <- sq("CTGA")
#'
#' # Using a method for character vector:
#' sq_2 <- as.sq("CTGA")
#'
#' # Checking that both objects are identical:
#' identical(sq_1, sq_2)
#'
#' @family output_functions
#' @export
as.sq <- function(x, ...)
  UseMethod("as.sq")

#' @rdname as.sq
#' @export
as.sq.default <- function(x, ...)
  import_sq(x, ...)

#' @rdname as.sq
#' @export
as.sq.character <- function(x, ...)
   sq(x, ...)

#' Convert sq object into character vector
#' 
#' @description Coerces sequences from an \code{\link[=sq-class]{sq}} object to
#' \code{\link{character}} vector of sequences.
#' 
#' @template x
#' @template NA_letter
#' @template three-dots
#' 
#' @return A \code{character} vector where each element represents the content
#' of respective sequence in input \code{sq} object.
#' 
#' @details
#' This method for \code{\link[=sq-class]{sq}} class allows converting sequences
#' from the \code{sq} object into a character vector of length equal to the
#' length of input. Each element of resulting vector is a separate sequence.
#' All attributes of the input sq are lost during the conversion to character
#' vector.
#'
#' @examples
#' # Creating an object to work on:
#' sq_dna <- sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", "CAGACCANNNATAG"),
#'              alphabet = "dna_ext")
#'
#' # Converting sq object into a character vector:
#' as.character(sq_dna)
#'
#' @family output_functions
#' @export
as.character.sq <- function(x, ...,
                            NA_letter = getOption("tidysq_NA_letter"))
  unpack(x, "STRING", NA_letter)

#' Convert sq object into matrix
#' 
#' @description Coerces sequences from a \code{\link[=sq-class]{sq}} object to
#' a \code{\link{matrix}}, in which rows correspond to sequences and columns to
#' positions.
#' 
#' @template x
#' @template three-dots
#' 
#' @return A \code{\link{matrix}} with number of rows the same as number of
#' sequences and number of columns corresponding to the length of the longest
#' sequence in the converted sq object.
#' 
#' @details
#' This method for class \code{sq} allows converting sequences from the
#' \code{sq} object into a matrix. Each row corresponds to the separate sequence
#' from the \code{sq} object, whereas each column indicates a single position
#' within a sequence. Dimensions of matrix are determined by the number of
#' sequences (rows) and the length of the longest sequence (columns). If length
#' of a sequence is smaller than the length of the longest sequence, the
#' remaining columns are filled with \code{NA}. All attributes of the input
#' \code{sq} are lost during the conversion to matrix.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("CGATAGACA", "TGACAAAAC", "GTGACCGTA"),
#'              alphabet = "dna_bsc")
#' sq_rna <- sq(c("CUGAAUGCAGUACCGUAAU", "AUGCCGUAAAUGCCAU", "CAGACCANNNAUAG"),
#'              alphabet = "rna_ext")
#'
#' # Sequences of the same lengths can be converted easily:
#' as.matrix(sq_dna)
#'
#' # Sequences that differ in length are filled with NA to the maximum length:
#' as.matrix(sq_rna)
#'
#' @family output_functions
#' @export
as.matrix.sq <- function(x, ...) {
  max_len <- max(get_sq_lengths(x))
  ret <- do.call(rbind, lapply(unpack(x, "STRINGS"), function(row) row[1:max_len]))
  ret[ret == getOption("tidysq_NA_letter")] <- NA
  ret
}

#' Check if object has specified type
#' 
#' @description Checks if object is an \code{\link[=sq-class]{sq}} object
#' without specifying type or if it is an \code{sq} object with specific type.
#'
#' @template x
#'
#' @return A \code{logical} value - \code{TRUE} if \code{x} has specified type,
#' \code{FALSE} otherwise.
#' 
#' @details
#' These functions are mostly simply calls to class checks. There are also
#' grouped checks, i.e. \code{is.sq_dna}, \code{is.sq_rna} and \code{is.sq_ami}.
#' These check for sq type regardless of if the type is basic or extended.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("GGCAT", "TATC-A", "TGA"), alphabet = "dna_bsc")
#' sq_rna <- sq(c("CGAUUACG", "UUCUAGA", "UUCA"), alphabet = "rna_bsc")
#' sq_ami <- sq(c("CVMPQGQQ", "AHLC--PPQ"), alphabet = "ami_ext")
#' sq_unt <- sq("BAHHAJJ&HAN&JD&", alphabet = "unt")
#' sq_atp <- sq(c("mALPVQAmAmA", "mAmAPQ"), alphabet = c("mA", LETTERS))
#'
#' # What is considered sq:
#' is.sq(sq_dna)
#' is.sq(sq_rna)
#' is.sq(sq_ami)
#' is.sq(sq_unt)
#' is.sq(sq_atp)
#'
#' # What is not:
#' is.sq(c(1,2,3))
#' is.sq(LETTERS)
#' is.sq(TRUE)
#' is.sq(NULL)
#'
#' # Checking for exact class:
#' is.sq_dna_bsc(sq_dna)
#' is.sq_dna_ext(sq_rna)
#' is.sq_rna_bsc(sq_ami)
#' is.sq_rna_ext(sq_rna)
#' is.sq_ami_bsc(sq_ami)
#' is.sq_ami_ext(sq_atp)
#' is.sq_atp(sq_atp)
#' is.sq_unt(sq_unt)
#'
#' # Checking for generalized type:
#' is.sq_dna(sq_atp)
#' is.sq_rna(sq_rna)
#' is.sq_ami(sq_ami)
#'
#' @family type_functions
#' @family util_functions
#' @export
is.sq <- function(x)
  test_class(x, "sq")

#' @rdname is.sq
#' @export
is.sq_dna_bsc <- function(x)
  test_class(x, "sq_dna_bsc")

#' @rdname is.sq
#' @export
is.sq_dna_ext <- function(x)
  test_class(x, "sq_dna_ext")

#' @rdname is.sq
#' @export
is.sq_dna <- function(x)
  is.sq_dna_bsc(x) || is.sq_dna_ext(x)

#' @rdname is.sq
#' @export
is.sq_rna_bsc <- function(x)
  test_class(x, "sq_rna_bsc")

#' @rdname is.sq
#' @export
is.sq_rna_ext <- function(x)
  test_class(x, "sq_rna_ext")

#' @rdname is.sq
#' @export
is.sq_rna <- function(x)
  is.sq_rna_bsc(x) || is.sq_rna_ext(x)

#' @rdname is.sq
#' @export
is.sq_ami_bsc <- function(x)
  test_class(x, "sq_ami_bsc")

#' @rdname is.sq
#' @export
is.sq_ami_ext <- function(x)
  test_class(x, "sq_ami_ext")

#' @rdname is.sq
#' @export
is.sq_ami <- function(x)
  is.sq_ami_bsc(x) || is.sq_ami_ext(x)

#' @rdname is.sq
#' @export
is.sq_unt <- function(x)
  test_class(x, "sq_unt")

#' @rdname is.sq
#' @export
is.sq_atp <- function(x)
  test_class(x, "sq_atp")

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
#' interpretation must be identical). If one of compared objects is a scalar,
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
#'  Multiple objects. For exact behavior, check Details section. First argument
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
#' class-less list. This behavior is effectively impossible and definitely
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
