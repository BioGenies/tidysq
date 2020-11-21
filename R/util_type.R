interpret_type <- function(name) {
  # TODO: improve; just improve
  switch(name, 
         dna_bsc = "dna_bsc",
         dna_ext = "dna_ext",
         rna_bsc = "rna_bsc",
         rna_ext = "rna_ext",
         ami_bsc = "ami_bsc",
         ami_ext = "ami_ext",
         unt = "unt")
}

type_as_class <- function(type)
  paste0("sq_", type)

#' Get type of a sq object
#' 
#' Function checks which type of sequences are contained in \code{\link{sq}} object.
#'  
#' @param x a \code{\link{sq}} object to be checked.
#'  
#' @return A \code{\link{character}} string, type of\code{\link{sq}} object - can be one of
#' "ami", "dna", "rna", "unt", "atp" or "enc".
#' 
#' @details This function returns type of sequence from \code{\link{sq}} object.
#' If the type of sequence is \strong{dna}, \strong{rna}, \strong{ami}, \strong{unt},
#' \strong{atp} or \strong{enc} function returns "dna", "rna", "ami", "unt", "atp" or
#' "enc" respectivetly.
#'  
#' @seealso \code{\link{sq}} \code{\link{sq}}
#' @export
get_sq_type <- function(x) {
  # TODO: a generic, maybe?
  assert_class(x, "sq")
  vec_ptype_abbr(x)
}
