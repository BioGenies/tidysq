#' sq: class for keeping biological sequences tidy
#' 
#' An object of class \strong{sq} represents a list of biological sequences. It is main
#' internal format of \strong{tidysq} package and most functions operate on it. 
#' The storage method is memory-optimized so that objects require as little memory
#' as possible (details below).
#' 
#' @section Construction/reading/import of sq objects:
#' There are multiple ways of obtaining \code{sq} objects:
#' \itemize{
#' \item constructing from a character vector with \code{\link{construct_sq}},
#' \item constructing from a character vector with \code{\link{as.sq}} method,
#' \item reading from the fasta file with \code{\link{read_fasta}},
#' \item importing from a format of other package like \code{ape} or
#' \code{Biostrings} with \code{\link{import_sq}}.
#' }
#'
#' \strong{Important note:} A manual assignment of a class \code{sq} to an object is
#' \strong{strongly discouraged} - due to the usage of low-level functions for
#' bit packing such assignment may lead to calling one of those functions during
#' operating on object or even printing it which can cause crash of R session and,
#' in consequence, loss of data.
#'
#' @section Export/writing of sq objects:
#' There are multiple ways of saving \code{sq} objects or converting them into
#' other formats:
#' \itemize{
#' \item converting into a character vector with \code{\link[sq:as.character.sq]{as.character}} method,
#' \item converting into a character matrix (or numeric matrix, in case of \strong{enc}) 
#' with \code{\link[sq:as.matrix.sq]{as.matrix}} method,
#' \item converting into a list of numerics with \code{\link{encsq_to_list}} (only
#' for \strong{enc} \code{sq}),
#' \item saving into the fasta file with \code{\link{write_fasta}},
#' \item exporting into a format of other package like \code{ape} or
#' \code{Biostrings} with \code{\link{export_sq}}.
#' }
#'
#' @section Types of sq:
#' This package is meant to handle both amino acids and nucleotides sequences 
#' thus there is need to differentiate \code{sq} objects that keep them. In 
#' addition, there are special types for handling non-standard sequence 
#' formats and encodings.
#' 
#' Each \strong{sq} object has exactly one of \strong{types}:
#' \itemize{
#' \item \strong{ami} - (\emph{amino acids}) represents a list of sequences of amino acids
#' (peptides or proteins),
#' \item \strong{nuc} - (\emph{nucleotides}) represents a list of nucleotide sequences (RNA or 
#' DNA),
#' \item \strong{unt} - (\emph{untyped}) represents a list of sequences that do not have 
#' specified type. They are mainly result of reading sequences from a file that 
#' contains some letters that are not in standard nucleotide or amino acid alphabets 
#' and user has not specified them explicitly. They should be converted to \strong{ami} 
#' or \strong{nuc} sequences (using functions like \code{\link{substitute_letters}} or 
#' \code{\link{typify}}).
#' \item \strong{atp} - (\emph{atypical}) represents sequences that have an alphabet
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
#' Additionally, there is a special subtype \strong{cln} (standing for \emph{clean}).
#' Only \strong{ami} and \strong{nuc} \code{sq} objects may have this subtype. It indicates
#' that sequences do not contain ambiguous letters (see "alphabets" section below).
#' 
#' \code{sq} object type is printed when using overloaded method 
#' \code{\link[sq:print.sq]{print}}. It can be also checked by using \code{\link{get_sq_type}}
#' @section Alphabet:
#' Each \code{sq} object have an \strong{alphabet} associated with it. Alphabet is
#' a set of possible \strong{letters} that can appear in sequences contained in object.
#' Alphabet is kept mostly as a character vector, where each element represents one
#' \strong{letter}.
#'
#' \code{sq} objects of type \strong{ami} or \strong{nuc} have fixed alphabets (depending
#' also if object has \strong{cln} subtype). In other words, if two \code{sq} objects have
#' exactly the same type - \strong{ami} or \strong{nuc} - and either both have or both don't
#' have \strong{cln} subtype, they are ensured to have the same alphabets.
#'
#' Here are listed alphabets for these types:
#' \itemize{
#' \item \strong{ami} \strong{cln} - ACDEFGHIKLMNPQRSTVWY-*
#' \item \strong{ami}, (not \strong{cln}) - ABCDEFGHIJKLMNOPQRSTUVWXYZ-*
#' \item \strong{nuc} \strong{cln} - ACGTU-
#' \item \strong{nuc}, (not \strong{cln}) - ACGTUWSMKRYBDHVN-
#' }
#'
#' To see details of these alphabets see \code{\link{aminoacids_df}} and 
#' \code{\link{nucleotides_df}}.
#'
#' Other types of \code{sq} objects are allowed to have different alphabets. Having an alphabet
#' exactly identical to one of those above does not automatically indicate that type of
#' sequence is one of those - e.g. there might be \strong{unt} \code{sq} that has an alphabet
#' identical to \strong{ami} \strong{cln} alphabet. To set the type, you should use the
#' \code{\link{typify}} function.
#'
#' The purpose of co-existence of \strong{unt} and \strong{atp} alphabets is the 
#' fact that although there is a standard for format of \emph{fasta} files, sometimes 
#' there are other types of symbols which do not fit the standard. Thanks to these types, 
#' tidysq can import files with customized alphabets can be imported. Moreover, an user 
#' may want to group amino acids with simillar properties (e.g. for machine learning) 
#' and replace the longer alphabet with symbols of groups. To check details, see 
#' \code{\link{read_fasta}}, \code{\link{construct_sq}} and \code{\link{substitute_letters}}.
#'
#' All of the types: \strong{ami}, \strong{nuc}, \strong{atp}, \strong{unt} have alphabets
#' that are character vectors while \strong{enc} objects have alphabets that are numeric
#' vectors. They are result of \code{\link{encode}} function, see its manual for details.
#'
#' \strong{Important note:} in \strong{atp} alphabets there is possibility of appearance of
#' letters that consists of more than one character - this functionality is provided in
#' order to handle situations like post-translational modifications, (e.g., using "mA" to 
#' indicating methylated alanine).
#'
#' \strong{Important note:} alphabets of \strong{atp} and \strong{unt} \code{sq} objects
#' are case sensitive. Thus, in their alphabets there can appear both lowercase and
#' uppercase letters simultaneously and they are treated as different characters. Alphabets of
#' \strong{nuc} and \strong{ami} objects are always uppercase and all functions converts
#' other parameters to uppercase when working with \strong{ami} or \strong{nuc} - e.g.
#' \link{\%has\%} operator converts lower letters to upper when searching for motifs in
#' \strong{ami} or \strong{nuc} object.
#' 
#' \strong{Important note:} maximum length of an alphabet is \strong{30 letters}. You are
#' not allowed to read fasta files or construct from character vectros that have more
#' than 30 distinct characters in sequences (with exception of reading or constructing
#' \strong{ami} or \strong{nuc} objects - during their construction lowercase letters
#' are automatically converted to uppercase).
#'
#' You can obtain an alphabet of the \code{sq} object using the \code{\link{get_sq_alphabet}}
#' function. You can
#' check which letters are invalid (are not represented in standard amino acids or nucleotides alphabet)
#' in each sequence of given \code{sq} object of type \strong{unt} or \strong{atp} by using
#' \code{\link{get_invalid_letters}}. You can substitute one letter with another using
#' \code{\link{substitute_letters}}.
#'
#' @section Missing/Not Available values:
#' There is a possibility of introducing \code{NA} values into sequences. \code{NA} value
#' does not represents gap (which are represented by \code{-}) or wildcard elements 
#' (\code{N} in he case of nucleotides and \code{X} in the case of amino acids), but is used 
#' as a representation of an empty position or invalid letters (not represented in nucleotide or 
#' amino acid alphabet).
#' \code{NA} does not belong to any alphabet (with exception of \strong{enc} objects where some of
#' letters might be \code{NA}, see \code{\link{encode}} for details). It is printed as
#' "!" and, thus, it is highly unrecommended to use "!" as special letter in \strong{atp}
#' sequences (but print character can be changed in options, see \code{\link{sq-options}}).
#'
#' \code{NA} might be introduced by:
#' \itemize{
#' \item reading fasta file with non standard letters in
#' \code{\link[sq:fast-mode]{fast mode}} with \code{\link{read_fasta}},
#' \item replacing a letter with \code{NA} value with \code{\link{substitute_letters}},
#' \item subsetting sequences out of their lengths with \code{\link{bite}}.
#' }
#'
#' An user can convert sequences that contain \code{NA} values into \code{NULL}
#' sequences with \code{\link{remove_na}}.
#' 
#' @section NULL (empty) sequences:
#' \code{NULL} sequence is a sequence of length 0.
#'
#' \code{NULL} sequences might be introduced by:
#' \itemize{
#' \item constructing \code{sq} object from character string of length zero with
#' \code{\link{construct_sq}},
#' \item using the \code{\link{clean}} function,
#' \item using the \code{\link{remove_na}} function,
#' \item subsetting \code{sq} object with \link[sq:sqextract]{extract operator}
#' }
#'
#' @section Storage format:
#' \code{sq} object is, in fact, \strong{list of raw vectors}. The fact that it is list
#' implies that an user can concatenate \code{sq} objects using \code{\link[sq:c.sq]{c}} method
#' and subset them using \code{\link[sq:sqextract]{extract operator}}. Alphabet is kept
#' as an attribute of the object. 
#' 
#' Raw vectors are the most efficient way of storage - each letter of sequece has asigned
#' an integer (its index in alphabet of \code{sq} object). Those integers in binary format
#' fit in less than 8 bits, but normally are stored on 16 bits. However, thanks to bit
#' packing it is possible to remove unused bits and store numbers more tightly. This 
#' operations result in a little time overhead in all operations, because most of them 
#' require unpacking and repacking sequences, but this cost is relatively low in comparision
#' to amount of saved memory.
#' 
#' For example - \strong{nuc} \strong{cln} alphabet consist of 6 values: ACGTU-. They are 
#' asigned numbers 1 to 6 respectively. Those numbers in binary format take form: \code{001},
#' \code{010}, \code{011}, \code{100}, \code{101}, \code{110}. Each of the letters can 
#' be coded with just 3 bits instead of 8 which is demanded by \code{char} - this allows
#' us to save more than 60\% of memory spent on storage of nucleotides sequences.
#' 
#' @section tibble compatibility:
#' \code{sq} objects are compatible with \code{tibble} class - that means you can have
#' \code{sq} object as a column of \code{tibble}. There are overloaded print methods, so
#' that it is printed in pretty format.
#'
#' @name sq
NULL

#' Construct sq object from character vector
#' 
#' @description This function allows the user to construct objects of 
#' \code{\link[sq:sq]{class sq}} from a character vector.
#' 
#' @param sq \code{\link{character}} vector 
#' @param type \code{\link{character}} string indicating type of \code{sq} object that
#' is going to be constructed; supported values are "ami" for amino acid sequences,
#' "nuc" for nucleotide sequences, "unt" for and \code{NULL} for type guessing (see details)
#' @param is_clean \code{\link{logical}} value indicating if sequences are clean,
#' or in other words - they don't contain ambiguous values; supported values are \code{TRUE} 
#' for clean sequences, \code{FALSE} for unclean sequences and \code{NULL} for auto detecting
#' (see details)
#' @param non_standard \code{\link{character}} vector indicating non-standard letters
#' contained in sequences. If \code{NULL}, sequences won't be searched for non-standard letters
#' of length more than one. Each element of this parameter should be at least two characters 
#' long
#' @return object of \code{\link[sq:sq]{class sq}} with appropriate type (one of: \strong{ami},
#' \strong{nuc}, \strong{unt}, \strong{atp}).
#' 
#' @details 
#' Function covers all possibilities of standard and non-standard types and alphabets.
#' You can check what 'type' and 'alphabet' exactly is in \code{\link{sq}} documentation.
#' Below there is a guide how function operates and how the program behaves according to the given 
#' arguments and the letters in the sequences.
#' 
#' \code{sq} parameter should be a character vector. Each element of this vector is a biological 
#' sequence. If this parameter has length 0, object of class \code{sq} with 0 sequences will be 
#' created (if not specified, it will have \strong{nuc} \strong{cln} type, which is a result of 
#' rules written below). If it contains sequences of length 0, \code{\link{NULL}} sequences
#' will be introduced (see \emph{NULL (empty) sequences} section in \code{\link[sq]{sq}}).
#' 
#' \strong{Important note:} in all below cases word 'letter' stands for an element of an alphabet.
#' Letter might consist of more than one character, for example "Ala" might be a single letter.
#' However, if you want to construct or read sequences with multi-character letters, you have 
#' to specify \code{non_standard} parameter. Details of letters, alphabets and types can be 
#' found in \code{\link[sq:sq]{sq class}} documentation.
#' 
#' @section Simple guide to constructing:
#' In most cases, all you need to do is just specifying \code{sq} parameter - type of sequences
#' will be guessed accordingly to rules described below. You need to pay attention, however, 
#' because for short sequences type may be guessed incorrectly - in this case you should
#' specify \code{type} and/or \code{is_clean}.
#' 
#' If your sequences contain non-standard letters, where each non-standard letter is one
#' character long, you also don't need to specify any parameter. Optionally, you can explicitly
#' do it by setting \code{type} to "atp".
#' 
#' If you want to construct sequences with multicharacter letters, you have to specify 
#' \code{non_standard} parameter, where you have to provide all non-standard letters longer
#' than one character.
#' 
#' In \code{\link[sq:fast-mode]{fast mode}} you have to specify both \code{type} and \code{is_clean} 
#' parameters. You cannot specify \code{non_standard} parameter in this mode. All letters
#' outside specified alphabet will be red as \code{\link[NA]{NA values}}.
#' 
#' @section Detailed guide to constructing:
#' Below there are listed all possibilities that can happen during constructing a \code{sq} object.
#' 
#' In normal mode (no \code{\link[sq:fast-mode]{fast mode}}):
#' \itemize{
#' \item If you don't specify any other parameter than \code{sq}, function will try to guess
#' sequence type and if it's clean (it will check in exactly this order):
#' \enumerate{
#' \item If it contains only ACGTU- letters, either lowercase or 
#' uppercase, type will be set to \strong{nuc} \strong{cln}.
#' \item If it contains any letters from 1. and additionally letters DEFHIKLMNPQRSVWY*, either
#' lowercase or uppercase, type will be set to \strong{ami} \strong{cln}. 
#' \item If it contains any letters from 1. and additionally letters WSMKRYBDHVN, either 
#' lowercase or uppercase, and does not contain any other letter, type will be set to 
#' \strong{nuc} without \strong{cln} subtype.
#' \item If it contains any letters from 1., 2., 3. and additionally letters JOUXZ, type will
#' be set to \strong{ami} without \strong{cln} subtype.
#' \item If it contains any letters that exceed all groups mentioned above, type will be set
#' to "unt".
#' }
#' \item If you specify \code{type} parameter as "ami" or "nuc" (and do not specify neither 
#' \code{is_clean} nor \code{non_standard}) type will be checked during construction: 
#' \enumerate{
#' \item If all letters in sequences fit the clean alphabet of given type, type will be set to 
#' given type with \strong{cln} subtype. 
#' \item If all letters in sequences fit the unclean alphabet of given type, type will be set to 
#' given type without \strong{cln} subtype.
#' \item If at least one of sequences contain at least one letter that is not element of unclean 
#' alphabet of provided type, an error will be thrown. 
#' }
#' \item If you specify both \code{type} and \code{is_clean}, function checks if letters
#' in sequences matches exactly specified alphabet (with capitalisation accuracy). If they do, 
#' type will be set to it. Otherwise, an error will be thrown.
#' \item If you specify \code{type} as "unt" and won't neither \code{is_clean} nor 
#' \code{non_standard}, type will be set to \strong{unt}. Letters won't be converted to uppercase, 
#' alphabet will consist of all letters found in sequences. 
#' \item If you do not sepcify neither \code{type} nor \code{is_clean} and specify 
#' \code{non_standard} parameter, which should be character vector where each element is at least 
#' two characters long, all strings as specified will be detected in sequences and treated as 
#' letters in constructed \strong{atp} \code{sq}.
#' \item All other combinations of parameters are incorrect.
#' }
#' 
#' In \code{\link[sq:fast-mode]{fast mode}} you have to specify \code{type} (it has to have either
#' "ami" or "nuc" value) and \code{is_clean} (\code{TRUE} or \code{FALSE}). You cannot specify
#' \code{non_standard}. All letters that aren't elements of destination alphabet (with a letter 
#' size accuracy) will be treated as \link[NA]{NA values}.
#' 
#' @section Handling with atp and unt sq and NA values:
#' You can convert letters into another using \code{\link{substitute_letters}} and then you
#' can use \code{\link{typify}} function to set of \code{sq} to \strong{ami} or \strong{nuc}.
#' If your sequences contain \code{NA} values, use \code{\link{remove_na}}
#' 
#' @examples 
#' # saving option:
#' previous_option <- getOption("tidysq_fast_mode")
#' 
#' #### constructing sq in normal mode:
#' ## setting an option:
#' options(tidysq_fast_mode = FALSE)
#' 
#' ## constructing sq without specyfiing type
#' # nuc cln sq
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG")) 
#' 
#' # ami cln sq
#' construct_sq(c("YQQPAVVM", "PQCFL"))
#' 
#' # ami cln sq can contain * - letter meaning end of translation:
#' construct_sq(c("MMDF*", "SYIHR*", "MGG*"))
#' 
#' # nuc sq
#' construct_sq(c("WHDHKYN", "GCYVCYU"))
#' 
#' # ami sq
#' construct_sq(c("XYOQWWKCNJLO"))
#' 
#' # unt sq - let's assume that one wants to mark some special element in sequnece with %
#' construct_sq(c("%%YAPLAA", "PLAA"))
#' 
#' ## constructing sq with type 
#' # all above examples will result in an identical if specified type as guessed
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "nuc") 
#' construct_sq(c("YQQPAVVM", "PQCFL"), "ami")
#' construct_sq(c("MMDF*", "SYIHR*", "MGG*"), "ami")
#' construct_sq(c("WHDHKYN", "GCYVCYU"), "nuc")
#' construct_sq(c("XYOQWWKCNJLO"), "ami")
#' construct_sq(c("%%YAPLAA", "PLAA"), "unt")
#' 
#' # you can force type other than guessed (if letters fit in the destination alphabet)
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "nuc", is_clean = FALSE)
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "ami")
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "ami", is_clean = FALSE)
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "unt")
#' 
#' ## constructing with non_standard specified
#' # in sequences below "mA" denotes methyled alanine - two characters are treated as single letter
#' construct_sq(c("LmAQYmASSR", "LmASMKLKFmAmA"), non_standard = "mA")
#' 
#' # reading sequences with three-letter names:
#' construct_sq(c("ProProGlyAlaMetAlaCys"), non_standard = c("Pro", "Gly", "Ala", "Met", "Cys"))
#' 
#' #### constructing in fast mode:
#' ## setting fast mode on
#' options(tidysq_fast_mode = TRUE)
#' 
#' # you cannot construct without specifying type
#' \dontrun{
#' construct_sq("CTGA")
#' }
#' construct_sq("CTGA", "nuc", TRUE)
#' 
#' # you cannot construct with specifying non_standard
#' \dontrun{
#' construct_sq("mAPQ", non_standard = "mA")
#' }
#' 
#' # letters other than in specified alphabet will be treated as NA:
#' construct_sq("NCTGCNA", "nuc", TRUE)
#' 
#' #### Other examples:
#' ## setting fast mode off again:
#' options(tidysq_fast_mode = FALSE)
#' 
#' # lowercase letters are converted to uppercase if detected type is "ami" or "nuc"
#' construct_sq(c("aTGc", "tcgTTA", "tt--AG"))
#' construct_sq(c("XYOqwwKCNJLo"))
#' 
#' # but not for "unt"
#' construct_sq(c("aAAaAA"), type = "unt")
#' 
#' # you can construct sq with length 0
#' construct_sq(character(0))
#' 
#' # and sq with empty sequences
#' construct_sq(c("AGTGGC", "", "CATGA", ""))
#' 
#' ## reseting an option
#' options(tidysq_fast_mode = previous_option)
#' 
#' @seealso \code{\link{sq}} \code{\link{read_fasta}} \code{\link{sq-options}} 
#' \code{\link{fast-mode}} \code{\link{substitute_letters}} \code{\link{remove_na}}
#' @exportClass sq
#' @export
construct_sq <- function(sq, type = NULL, is_clean = NULL, non_standard = NULL) {
  .check_character(sq, "'sq'", allow_zero_len = TRUE)
  if (.is_fast_mode()) {
    .nc_construct_sq(sq, type, is_clean)
  } else {
    .check_type_or_nonst_alph(type, is_clean, non_standard)
    if (!is.null(non_standard)) {
      .nonst_construct_sq(sq, non_standard)
    } else {
      .check_logical(is_clean, "'is_clean'", allow_null = TRUE, single_elem = TRUE)
      if (is.null(type)) {
        type_clean <- .guess_sq_type_subtype(sq)
        type <- type_clean[["type"]]
        if (is.null(is_clean) && type != "unt") is_clean <- type_clean[["is_clean"]]
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
      } else {
        res_ind <- apply(binded, 1, function(row) row[1]:row[2])
        res_ind <- if (is.list(res_ind)) unlist(res_ind) else as.integer(res_ind)
      }
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