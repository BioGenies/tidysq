#' Remove sequences containing NA values
#' 
#' Removes sequences containing ambiguous elements or removes \code{\link[=sq]{NA values}} 
#' from sequences in a \code{\link{sq}} object.
#' 
#' @inheritParams clean
#'  
#' @return A \code{\link{sq}} object with the same type as input type. Sequences not containing
#' any \code{\link[=sq]{NA}} values are left untouched.
#' 
#' @details This function allows removal of sequences containing \code{\link[=sq]{NA}} values.
#' By default, whole sequences containing ambiguous elements are removed 
#' and \code{\link[=sq]{NULL}} (empty) sequences are introduced in their place. If 
#' \code{only_elements = TRUE} then only \code{\link[=sq]{NA}} values are removed 
#' from sequences in \code{sq} object. \code{\link[=sq]{NULL}} values (empty sequences) 
#' can be identified using \code{\link{is_null_sq}} function. 
#' 
#' \code{NA} may be introduced as a result of using functions like 
#' \code{\link{substitute_letters}} or \code{\link{bite}}. They also appear in sequences if
#' you are reading file using \code{\link{read_fasta}} or constructing \code{sq} object from
#' \code{\link{character}} vector with \code{\link{construct_sq}} in 
#' \code{\link[=fast-mode]{fast mode}} and there are letters in file or in strings other than
#' specified.
#'
#' @examples 
#' # Creating objects to work on:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA", "ACTNNAGCN"), type = "dna")
#' 
#' # Substituting some letters with NA
#' sq_ami_sub <- substitute_letters(sq_ami, c(E = NA, R = NA))
#' sq_dna_sub <- substitute_letters(sq_dna, c(N = NA))
#' 
#' # Biting sequences out of range
#' sq_bitten <- bite(sq_ami, 1:15)
#' 
#' # Printing them
#' sq_ami_sub
#' sq_dna_sub
#' 
#' # Removing sequences containing NA
#' remove_na(sq_ami_sub)
#' remove_na(sq_dna_sub)
#' remove_na(sq_bitten)
#' 
#' # Removing only NA elements
#' remove_na(sq_ami_sub, only_elements = TRUE)
#' remove_na(sq_dna_sub, TRUE)
#' remove_na(sq_bitten, TRUE)
#' 
#' @seealso \code{\link{sq}} \code{\link{is_null_sq}} \code{\link{substitute_letters}}
#' \code{\link{bite}}
#' @export
remove_na <- function(sq, only_elements = FALSE) {
  .validate_sq(sq)
  .check_logical(only_elements, "'only_elements'", single_elem = TRUE)
  
  alph <- alphabet(sq)
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  
  if (only_elements) {
    ret <- .apply_sq(sq, "int", "int", function(s) {
      s[s != na_val]
    }) 
  } else {
    ret <- lapply(sq, function(s) {
      st <- C_unpack_ints(s, alph_size)
      if (any(st == na_val)) structure(raw(), original_length = 0) else s
    })
  }
  vec_restore(ret, sq)
}