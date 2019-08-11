#' sq: class for keeping sequences in a tidy format
#' 
#' Object of class \strong{sq} represents list of biological sequences. It is main
#' internal format of \strong{tidysq} package and most functions operate on them. 
#' Storage method is memory-optimized so that objects require as little memory
#' as possible (details below).
#' 
#' @section Construction/reading/import of sq objects:
#' There are multiple ways of obtaining \code{sq} objects:
#' \itemize{
#' \item constructing from character vector with \code{\link{construct_sq}},
#' \item constructing from character vector with \code{\link{as.sq}} method,
#' \item reading from fasta file with \code{\link{read_fasta}},
#' \item exporting from format of other package like \code{ape} or
#' \code{Biostrings} with \code{\link{import_sq}}.
#' }
#'
#' \strong{Important note:} Manual assignment of a class \code{sq} to an object is
#' \strong{strongly discouraged} - due to usage of low-level functions for
#' bit packing such assignment may lead to calling one of those functions during
#' operating on object or even printing it which can cause crash of R session and,
#' in consequence, loss of data.
#'
#' @section Types of sq:
#' This package is meant to handle both amino acids and nucleotides sequences 
#' thus there is need to differentiate \code{sq} objects that keep them. In 
#' addition, there are special types for handling non-standard sequence 
#' formats and encodings.
#' 
#' Each \strong{sq} object has exactly one of \strong{types}:
#' \itemize{
#' \item \strong{ami} - (\emph{amino acids}) represents list of amino acids residues 
#' sequences
#' \item \strong{nuc} - (\emph{nucleotides}) represents list of nucleotides sequences
#' \item \strong{unt} - (\emph{untyped}) represents list of sequences that don't have 
#' specified type. They are mainly result of reading sequences from a file that 
#' contains some letters that aren't in standard alphabets and user hasn't 
#' specified them explicitly. They should be converted to \strong{ami} or \strong{nuc} 
#' sequences before beggining of work with them (using functions like 
#' \code{\link{substitute_letters}} or \code{\link{typify}}).
#' \item \strong{atp} - (\emph{atypical}) represents sequences that have alphabet
#' different from standard \strong{ami} or \strong{nuc} alphabets - similarly to 
#' \strong{unt}, but user has explicitly informed about it. They are
#' result of constructing sequences or reading from file with specifying 
#' \code{non_standard} parameter (for details see \code{\link{read_fasta}}
#' and \code{\link{construct_sq}}). They are also result of using function
#' \code{\link{substitute_letters}} - user can use this to for example simplify 
#' alphabet and replace a few letters with one.
#' \item \strong{enc} - (\emph{encoded}) represents list of sequences that have been
#' encoded with function \code{\link{encode}} where each letter is assigned with 
#' numeric value.
#' }
#' 
#' Additionally, there is special subtype \strong{cln} (standing for \emph{clean}).
#' Only \strong{ami} and \strong{nuc} \code{sq} objects may have this subtype. It indicates
#' that sequences don't contain ambiguous letters (see "alphabets" section below).
#' 
#' \code{sq} object type is printed when using overloaded method 
#' \code{\link[print.sq]{print}}. It can be also checked by using \code{\link{get_sq_type}}
#' @section Alphabet:
#' Each \code{sq} object have an \strong{alphabet} associated with it. Alphabet is
#' a set of possible \strong{letters} that can appear in sequences contained in object.
#' Alphabet is kept mostly as a character vector, where each element represents one
#' \strong{letter}.
#'
#' \code{sq} obejcts of type \strong{ami} or \strong{nuc} have fixed alphabets (depending
#' also if object has \strong{cln} subtype). In other words, if two \code{sq} objects have
#' exactly the same type - \strong{ami} or \strong{nuc} - and either both have or both don't
#' have \strong{cln} subtype, they are ensured to have the same alphabets.
#'
#' Here are listed alphabets for these types:
#' \itemize{
#' \item \strong{ami} (not \strong{cln}) - ACDEFGHIKLMNPQRSTVWY-*
#' \item \strong{ami}, \strong{cln} - ABCDEFGHIJKLMNOPQRSTUVWXYZ-*
#' \item \strong{nuc} (not \strong{cln}) - ACGTU-
#' \item \strong{nuc}, \strong{cln} - ACGTUWSMKRYBDHVN-
#' }
#'
#' To see details of these alphabets or find out how to obtain them through code, see
#' \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}}.
#'
#' Other types of \code{sq} objects are allowed to have other alphabets. Having alphabet
#' exactly identical to one of those above don't automatically indicates that type of
#' sequence is one of those - e.g. there might be \strong{unt} \code{sq} that has alphabet
#' identical to \strong{ami} \strong{cln} alphabet. To set the type, you should use
#' \code{\link{typify}} function.
#'
#' Purpose of existing \strong{unt} and \strong{atp} alphabets is the fact that although
#' there is a standard for format of \emph{fasta} files, sometimes there are another types
#' of symbols which doesn't fit the standard. Thanks to these types, files with customized
#' alphabets can be red. What's more, user may want to group amino acids with simillar
#' properties (e.g. for machine learning) and replace them with symbols of groups.
#' To check details, see \code{\link{read_fasta}}, \code{\link{construct_sq}} and
#' \code{\link{substitute_letters}}.
#'
#' All of types: \strong{ami}, \strong{nuc}, \strong{atp}, \strong{unt} have alphabets
#' that are character vectors while \strong{enc} objects have alphabets that are numeric
#' vectors. They are result of \code{\link{encode}} function, see its manual for details.
#'
#' \strong{Important note:} in \strong{atp} alphabets there is possibility of appearance of
#' letters that consists of more than one character - this functionality is provided in
#' order to handle situations like, for example, using "mA" to indicating methylated
#' alanine.
#'
#' \strong{Important note:} alphabets of \strong{atp} and \strong{unt} \code{sq} objects
#' are case sensitive. That means, in their alphabets there can appear both lowercase and
#' uppercase letters simultaneously and they are treated as different. Alphabets of
#' \strong{nuc} and \strong{ami} objects are uppercase and all functions converts
#' other parameters to uppercase when working with \strong{ami} or \strong{nuc} - e.g.
#' \link{`\%has\%`} operator converts lower letters to upper when searching for motifs in
#' \strong{ami} or \strong{nuc} object.
#'
#' You can obtain alphabet of \code{sq} object using \code{\link{get_sq_alphabet}}. You can
#' check which letters are invalid (aren't in standard amino acids or nucleotides alphabet)
#' in each sequence of given \code{sq} object of type \strong{unt} or \strong{atp} by using
#' \code{\link{get_invalid_letters}}. You can substitute one letter with another using
#' \code{\link{substitute_letters}}
#'
#' @section Missing/Not Available values:
#' There is a possibility of introducing \code{NA} values into sequences. \code{NA} value
#' represents element of sequence that is unknown, missing or invalid for some reason.
#' It isn't letter of any alphabet (with exception of \strong{enc} objects where some of
#' letters might be \code{NA}, see \code{\link{encode}} for details). It is printed as
#' "!" and, thus, it is highly unrecommended to use "!" as special letter in \strong{atp}
#' sequences (but print character can be changed in options, see \code{\link{sq-options}}).
#'
#' Missing values might be caused by:
#' \itemize{
#' \item reading fasta file with non standard letters in
#' \code{\link[no-check-mode]{no-check mode}} with \code{\link{read_fasta}},
#' \item replacing a letter with \code{NA} value with \code{\link{substitute_letters}},
#' \item subsetting sequences out of their lengths with \code{\link{bite}}.
#' }
#'
#' User shouldn't operate on sequences with \code{NA} values and should get rid of them
#' if possible. You can convert sequences that contain \code{NA} values into \code{NULL}
#' sequences with \code{\link{remove_na}}.
#' 
#' @section NULL (empty) sequences:
#' There is a possibility of appearance of \code{NULL} sequences in \code{sq} object.
#' \code{NULL} sequence is a sequence of length 0.
#'
#' \code{NULL} sequences might be caused by:
#' \itemize{
#' \item constructing \code{sq} object from character string of length zero with
#' \code{\link{construct_sq}},
#' \item using \code{\link{clean}} function,
#' \item using \code{\link{remove_na}} function,
#' \item subsetting \code{sq} object with \link[sq]{extract operator}
#' }
#'
#' @name sq
NULL


#' @exportClass sq
#' @export
construct_sq <- function(sq, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_character(sq, "'sq'", allow_zero_len = TRUE)
  if (getOption("tidysq_no_check_mode") == TRUE) {
    .nc_construct_sq(sq, type, is_clean)
  } else {
    .check_type_or_nonst_alph(type, is_clean, non_standard)
    if (!is.null(non_standard)) {
      .nonst_construct_sq(sq, non_standard)
    } else {
      .check_logical(is_clean, "'is_clean'", allow_null = TRUE, single_elem = TRUE)
      if (is.null(type)) {
        type <- .guess_sq_type(sq)
      }
      .check_type(type, allow_unt = TRUE)
      switch(type,
             ami = .construct_amisq(sq, is_clean),
             nuc = .construct_nucsq(sq, is_clean),
             unt = .construct_untsq(sq))
    }
  }
}

.nc_construct_sq <- function(sq, type, is_clean) {
  .check_type(type)
  .check_logical(is_clean, single_elem = TRUE)
  
  sq <- .nc_bitify_sq(sq, type, is_clean)
  sq <- .set_class(sq, type, is_clean)
  .set_alph(sq, .get_standard_alph(type, is_clean))
}


#' @importFrom stringi stri_sub
#' @importFrom stringi stri_locate_all_regex
.nonst_construct_sq <- function(sq, non_standard) {
  .check_character(non_standard, "'non_standard'")
  .check_nchar(non_standard, "'non_standard'", minimal_nchar = 2)
  
  sq <- lapply(sq, function(s) {
    pos <- stri_locate_all_regex(s, non_standard)
    binded <- apply(na.omit(do.call(rbind, pos)), 2, sort)
    if (length(binded) == 0) {
      strsplit(s, "")[[1]]
    } else {
      if (length(binded) == 2) {
        binded <- t(as.matrix(binded))
        res_ind <- binded[1]:binded[2]
      } else res_ind <- as.integer(binded, 1, function(row) row[1]:row[2])
      sin_ind <- setdiff(1:nchar(s), res_ind)
      ind <- .merge_ind(sin_ind, binded[,1])
      n <- nrow(binded) + length(sin_ind)
      begs <- integer(n)
      ends <- integer(n)
      begs[(1:n)[ind]] <- sin_ind
      ends[(1:n)[ind]] <- sin_ind
      begs[(1:n)[!ind]] <- binded[,1]
      ends[(1:n)[!ind]] <- binded[,2]
      stri_sub(s, begs, ends)
    }
  })
  
  alph <- unique(unlist(sq))
  .check_alph_length(alph)
  sq <- .bitify_sq(sq, alph)
  sq <- .set_alph(sq, alph)
  .set_class(sq, "atp")
}

#' @exportClass amisq
.construct_amisq <- function(sq, is_clean) {
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  .check_real_alph_clean(real_alph, "ami", is_clean)
  if (is.null(is_clean)) {
    is_clean <- .guess_ami_is_clean(real_alph)
  }
  sq <- .nc_bitify_sq(sq, "ami", is_clean)
  sq <- .set_alph(sq, .get_standard_alph("ami", is_clean))
  .set_class(sq, "ami", is_clean)
}

#' @exportClass nucsq
.construct_nucsq <- function(sq, is_clean) {
  sq <- toupper(sq)
  real_alph <- .get_real_alph(sq)
  .check_real_alph_clean(real_alph, "nuc", is_clean)
  if (is.null(is_clean)) {
    is_clean <- .guess_nuc_is_clean(real_alph)
  }
  sq <- .nc_bitify_sq(sq, "nuc", is_clean)
  sq <- .set_alph(sq, .get_standard_alph("nuc", is_clean))
  .set_class(sq, "nuc", is_clean)
}

#' @exportClass untsq
.construct_untsq <- function(sq) {
  alph <- .get_real_alph(sq)
  .check_alph_length(alph)
  
  sq <- .bitify_sq(sq, alph)
  sq <- .set_alph(sq, alph)
  .set_class(sq, "unt", FALSE)
}

validate_sq <- function(object, type = NULL) {
  argname <- deparse(substitute(object))
  if (!"sq" %in% class(object))
    stop(argname, " doesn't inherit class 'sq'", call. = FALSE)
  sqtype <- .get_sq_subclass(object)
  if (length(sqtype) != 1) 
    stop("'object' should have exactly one of classes: 'amisq', 'nucsq', 'untsq', 'encsq', 'atpsq'", call. = FALSE)
  if (is.null(attr(object, "alphabet"))) 
    stop("'object' doesn't 'alphabet' attribute", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.character(alph) &&
      !is.numeric(alph))
    stop("attribute 'alphabet' is neither a character nor a numeric vector", call. = FALSE)
  if (!is.list(object))
    stop("'object' isn't a list", call. = FALSE)
  if (!all(sapply(object, is.raw))) 
    stop("'object' isn't a list of raw vectors", call. = FALSE)
  if (!is.null(type)) {
    switch(type,
           ami = .validate_amisq(object),
           nuc = .validate_nucsq(object),
           unt = .validate_untsq(object),
           atp = .validate_atpsq(object),
           enc = .validate_encsq(object),
           invisible(object)
    )
  }
}

.validate_nucsq <- function(object) {
  if (!"nucsq" %in% class(object))
    stop("'object' doesn't inherit class 'nucsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!"clnsq" %in% class(object)) {
    if (!identical(alph, .get_standard_alph("nuc", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard nucleotides alphabet", call. = FALSE)
  } else if (!identical(alph, .get_standard_alph("nuc", TRUE)))
      stop("attribute 'alphabet' isn't identical to cleaned nucleotides alphabet", call. = FALSE)
  
  invisible(object)
}

.validate_amisq <- function(object) {
  if (!"amisq" %in% class(object)) 
    stop("'object' doesn't inherit class 'amisq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!"clnsq" %in% class(object)) {
    if (!identical(alph, .get_standard_alph("ami", FALSE)))
      stop("attribute 'alphabet' isn't identical to standard aminoacids alphabet", call. = FALSE)
  } else if (!identical(alph, .get_standard_alph("ami", TRUE))) 
      stop("attribute 'alphabet' isn't identical to cleaned aminoacids alphabet", call. = FALSE)
  
  
  invisible(object)
}

.validate_untsq <- function(object) {
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector", call. = FALSE)
  if (!"untsq" %in% class(object))
    stop("'object' doesn't inherit class 'untsq'", call. = FALSE)
  
  invisible(object)
}

.validate_atpsq <- function(object) {
  if (!"atpsq" %in% class(object))
    stop("'object' doesn't inherit class 'atpsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.character(alph))
    stop("attribute 'alphabet' isn't a character vector", call. = FALSE)
  
  invisible(object)
}

.validate_encsq <- function(object) {
  if (!"encsq" %in% class(object))
    stop("'object' doesn't inherit class 'encsq'", call. = FALSE)
  alph <- .get_alph(object)
  if (!is.numeric(alph))
    stop("attribute 'alphabet' isn't a numeric vector", call. = FALSE)
  
  invisible(object)
}