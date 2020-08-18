#' Extract parts of a sq object
#' 
#' @name sqextract
#' @aliases sq-extract
#' @description Operator to extract subsets of sq objects.
#' 
#' @param x sq object from which to extract element(s)
#' @param i,j,... indices specifying elements to extract. They may be 
#' \code{\link{numeric}}, \code{\link{character}} or \code{\link{logical}} vectors or empty. 
#' This function follows \code{\link[vctrs]{vctrs}} conventions regarding argument
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
NULL

#' Concatenate sq objects
#' 
#' @name sqconcatenate
#' @aliases sq-concatenate
#' @description Merges multiple \code{\link{sq}} objects of the same type into one larger.
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
NULL

#' Print sq object
#' 
#' @name sqprint
#' @aliases sq-print
#' @description Prints input \code{\link{sq}} object in a human-friendly form.  
#' 
#' @details \code{Print} method is used by default in each case of calling the 
#' \code{\link{sq}} object with default parameters. 
#' Only by explicit calling the \code{print} method parameters can be changed. 
#'  
#' \code{Print} checks if the input \code{\link{sq}} object is cleaned and includes 
#' this information alongside with type in the printed message. On the right side of 
#' the sequence, in angle brackets, the length of each sequence is printed (e.q. "<9>").
#' 
#' If the \code{max_sequences} parameter is supplied, the desired number of sequences 
#' is printed and this information is included in a message (e.q. "printed 1 out of 3"). 
#' Only \code{max_sequences} value smaller than the number of sequences in object 
#' affects the function. The default value indicating how many sequences should 
#' be printed is 10, but it can be changed in \code{\link[=tidysq-options]{package options}}. 
#' 
#' Default value of \code{use_color} parameter is \code{TRUE} - sequences are printed
#' in green, while empty sequences, NA character and dots in gray. If this option is disabled, 
#' all sequences are in default color of console.
#' 
#' The \code{letters_sep} parameter indicates how the letters should be separated 
#' (they are not by default). Any character string can be supplied but 
#' \code{\link{NA_character_}}.
#' 
#' If sequences are too long, only leading characters are printed (as many as possible
#' in single line) and following dots indicating that sequence is truncated.
#' 
#' If sequences contain \code{\link{NA}} (‘Not Available’ / Missing Values) values, they 
#' are printed as "!" character, but it can be changed in 
#' \code{\link[=tidysq-options]{package options}}.
#' 
#' This is overloaded function from base package. It is selected when \code{\link{sq}} 
#' object is used as a parameter for print function. To see the generic function 
#' page, check \code{\link[base:print]{here}}.
#' 
#' @param x \code{\link{sq}} object
#' @param max_sequences \code{numeric} value indicating how many sequences 
#' should be printed
#' @param use_color \code{logical} value indicating if sequences should 
#' be colored
#' @param letters_sep \code{character} value indicating how the letters 
#' should be separated
#' @param ... 	further arguments passed to or from other methods. 
#' Unused.
#' 
#' @examples
#' 
#' # Creating sq objects using construct_sq:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA"), type = "dna")
#' sq_unt <- construct_sq(c("ATGCAGGA!", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#' 
#' # Printing without explicit function calling with default parameters:
#' sq_ami
#' sq_dna
#' sq_unt
#' 
#' # Printing with explicit function calling and specific parameters:
#' print(sq_ami)
#' print(sq_dna, max_sequences = 1, use_color = FALSE)
#' print(sq_unt, letters_sep = ":")
#' 
#' # Printing of the cleaned object:
#' clean(sq_dna)
#' print(clean(sq_dna), letters_sep = "-", use_color = FALSE)
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{tidysq-options}}
NULL
