#' sq: class for keeping biological sequences tidy
#' 
#' An object of class \strong{sq} represents a list of biological sequences. It is the 
#' main internal format of the \pkg{tidysq} package and most functions operate on it. 
#' The storage method is memory-optimized so that objects require as little memory
#' as possible (details below).
#' 
#' @section Construction/reading/import of sq objects:
#' There are multiple ways of obtaining \code{sq} objects:
#' \itemize{
#' \item constructing from a \code{\link{character}} vector with \code{\link{construct_sq}},
#' \item constructing from another object with \code{\link{as.sq}} method,
#' \item reading from the fasta file with \code{\link{read_fasta}},
#' \item importing from a format of other package like \pkg{ape} or
#' \pkg{Biostrings} with \code{\link{import_sq}}.
#' }
#'
#' \strong{Important note:} A manual assignment of a class \code{sq} to an object is
#' \strong{strongly discouraged} - due to the usage of low-level functions for
#' bit packing such assignment may lead to calling one of those functions during
#' operating on object or even printing it which can cause a crash of R session and,
#' in consequence, loss of data.
#'
#' @section Export/writing of sq objects:
#' There are multiple ways of saving \code{sq} objects or converting them into
#' other formats:
#' \itemize{
#' \item converting into a character vector with \code{\link[=as.character.sq]{as.character}} 
#' method,
#' \item converting into a character matrix (or numeric matrix, in case of \strong{enc}) 
#' with \code{\link[=as.matrix.sq]{as.matrix}} method,
#' \item converting into a list of numerics with \code{\link{encsq_to_list}} (only
#' for \strong{enc} \code{sq}),
#' \item saving as fasta file with \code{\link{write_fasta}},
#' \item exporting into a format of other package like \code{ape} or
#' \code{Biostrings} with \code{\link{export_sq}}.
#' }
#'
#' @section Types of sq:
#' This package is meant to handle both amino acid, DNA and RNA sequences, 
#' thus there is need to differentiate \code{sq} objects that keep them. In 
#' addition, there are special types for handling non-standard sequence 
#' formats and encodings.
#' 
#' Each \strong{sq} object has exactly one of \strong{types}:
#' \itemize{
#' \item \strong{ami} - (\emph{amino acids}) represents a list of sequences of amino acids
#' (peptides or proteins),
#' \item \strong{dna} - (\emph{DNA}) represents a list of DNA sequences,
#' \item \strong{rna} - (\emph{RNA}) represents a list of RNA sequences (together with
#' \strong{dna} above often collectively called "nucleotide sequences"),
#' \item \strong{unt} - (\emph{untyped}) represents a list of sequences that do not have 
#' specified type. They are mainly result of reading sequences from a file that 
#' contains some letters that are not in standard nucleotide or amino acid alphabets 
#' and user has not specified them explicitly. They should be converted to \strong{ami} ,
#' \strong{dna} or \strong{rna} sequences (using functions like \code{\link{substitute_letters}}
#' or \code{\link{typify}}).
#' \item \strong{atp} - (\emph{atypical}) represents sequences that have an alphabet
#' different from standard \strong{ami}, \strong{dna} or \strong{rna} alphabets - similarly
#' to \strong{unt}, but user has explicitly informed about it. They are
#' result of constructing sequences or reading from file with specifying 
#' \code{non_standard} parameter (for details see \code{\link{read_fasta}}
#' and \code{\link{construct_sq}}). They are also result of using function
#' \code{\link{substitute_letters}} - users can use this for example to simplify 
#' an alphabet and to replace a few letters by one.
#' \item \strong{enc} - (\emph{encoded}) represents list of sequences that have been
#' encoded with function \code{\link{encode}} where each letter is assigned with 
#' numeric value.
#' }
#' 
#' Additionally, there is a special subtype \strong{cln} (standing for \emph{clean}).
#' Only \strong{ami}, \strong{dna} or \strong{rna} \code{sq} objects may have this subtype.
#' It indicates that sequences do not contain ambiguous letters (see "alphabets" section below).
#' 
#' \code{sq} object type is printed when using overloaded method 
#' \code{\link[=sq-print]{print}}. It can be also checked by using \code{\link{get_sq_type}}
#' @section Alphabet:
#' Each \code{sq} object have an \strong{alphabet} associated with it. Alphabet is
#' a set of possible \strong{letters} that can appear in sequences contained in object.
#' Alphabet is kept mostly as a character vector, where each element represents one
#' \strong{letter}.
#'
#' \code{sq} objects of type \strong{ami}, \strong{dna} or \strong{rna} have fixed alphabets
#' (depending also if object has \strong{cln} subtype). In other words, if two \code{sq}
#' objects have exactly the same type - \strong{ami}, \strong{dna} or \strong{rna} - and
#' either both have or both do not have a \strong{cln} subtype, they are ensured to have
#' the same alphabets.
#'
#' Here are listed alphabets for these types:
#' \itemize{
#' \item \strong{ami} \strong{cln} - ACDEFGHIKLMNPQRSTVWY-*
#' \item \strong{ami}, (not \strong{cln}) - ABCDEFGHIJKLMNOPQRSTUVWXYZ-*
#' \item \strong{dna} \strong{cln} - ACGT-
#' \item \strong{dna}, (not \strong{cln}) - ACGTWSMKRYBDHVN-
#' \item \strong{rna} \strong{cln} - ACGU-
#' \item \strong{rna}, (not \strong{cln}) - ACGUWSMKRYBDHVN-
#' }
#'
#' To see details of these alphabets see \code{\link{aminoacids_df}} and 
#' \code{\link{nucleotides_df}}.
#'
#' Other types of \code{sq} objects are allowed to have different alphabets. Having an alphabet
#' exactly identical to one of those above does not automatically indicate that the type of
#' the sequence is one of those - e. g., there might be \strong{unt} \code{sq} that has an alphabet
#' identical to \strong{ami} \strong{cln} alphabet. To set the type, one should use the
#' \code{\link{typify}} function.
#'
#' The purpose of co-existence of \strong{unt} and \strong{atp} alphabets is the 
#' fact that although there is a standard for format of \emph{fasta} files, sometimes 
#' there are other types of symbols, which do not match the standard. Thanks to these types, 
#' tidysq can import files with customized alphabets can be imported. Moreover, a user 
#' may want to group amino acids with similar properties (e. g., for machine learning) 
#' and replace the longer alphabet with symbols of groups. To check details, see 
#' \code{\link{read_fasta}}, \code{\link{construct_sq}} and \code{\link{substitute_letters}}.
#'
#' All of the types: \strong{ami}, \strong{dna}, \strong{rna}, \strong{atp}, \strong{unt}
#' have alphabets that are character vectors while \strong{enc} objects have alphabets that
#' are numeric vectors. They are result of \code{\link{encode}} function, see its manual
#' for details.
#'
#' \strong{Important note:} in \strong{atp} alphabets there is possibility of appearance of
#' letters that consists of more than one character - this functionality is provided in
#' order to handle situations like post-translational modifications, (e.g., using "mA" to 
#' indicating methylated alanine).
#'
#' \strong{Important note:} alphabets of \strong{atp} and \strong{unt} \code{sq} objects
#' are case sensitive. Thus, in their alphabets both lowercase and uppercase letters can
#' appear simultaneously and they are treated as different characters. Alphabets of
#' \strong{dna}, \strong{rna} and \strong{ami} objects are always uppercase and all functions converts
#' other parameters to uppercase when working with \strong{dna}, \strong{rna} or \strong{ami} - e.g.,
#' \code{\link{\%has\%}} operator converts lower letters to upper when searching for motifs in
#' \strong{dna}, \strong{rna} or \strong{ami} object.
#' 
#' \strong{Important note:} maximum length of an alphabet is \strong{30 letters}. You are
#' not allowed to read fasta files or construct from character vectors that have more
#' than 30 distinct characters in sequences (with exception of reading or constructing
#' \strong{ami}, \strong{dna} or \strong{rna} objects - during their construction lowercase
#' letters are automatically converted to uppercase).
#'
#' You can obtain an alphabet of the \code{sq} object using the \code{\link{get_sq_alphabet}}
#' function. You can check which letters are invalid (are not represented in standard
#' amino acid or nucleotide alphabet) in each sequence of given \code{sq} object of type
#' \strong{unt} or \strong{atp} by using \code{\link{get_invalid_letters}}. You can
#' substitute one letter with another using \code{\link{substitute_letters}}.
#'
#' @section Missing/Not Available values:
#' There is a possibility of introducing \code{\link{NA}} values into sequences. \code{NA} value
#' does not represents gap (which are represented by \code{'-'}) or wildcard elements 
#' (\code{N} in the case of nucleotides and \code{'X'} in the case of amino acids), but is used 
#' as a representation of an empty position or invalid letters (not represented in nucleotide or 
#' amino acid alphabet).
#' \code{NA} does not belong to any alphabet (with exception of \strong{enc} objects where some of
#' letters might be \code{NA}, see \code{\link{encode}} for details). It is printed as
#' '!' and, thus, it is highly unrecommended to use '!' as special letter in \strong{atp}
#' sequences (but print character can be changed in options, see \code{\link{tidysq-options}}).
#'
#' \code{NA} might be introduced by:
#' \itemize{
#' \item reading fasta file with non-standard letters in
#' \code{\link[=fast-mode]{fast mode}} with \code{\link{read_fasta}},
#' \item replacing a letter with \code{NA} value with \code{\link{substitute_letters}},
#' \item subsetting sequences beyond their lengths with \code{\link{bite}}.
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
#' \item subsetting \code{sq} object with \link[=sq-extract]{extract operator}
#' }
#'
#' @section Storage format:
#' \code{sq} object is, in fact, \strong{list of raw vectors}. The fact that it is list
#' implies that an user can concatenate \code{sq} objects using \code{\link[=sq-concatenate]{c}}
#' method and subset them using \code{\link[=sq-extract]{extract operator}}. Alphabet is kept
#' as an attribute of the object.
#' 
#' Raw vectors are the most efficient way of storage - each letter of the sequence has assigned
#' an integer (its index in alphabet of \code{sq} object). Those integers in binary format
#' fit in less than 8 bits, but normally are stored on 16 bits. However, thanks to bit
#' packing it is possible to remove unused bits and store numbers more tightly. This 
#' operations result in a little time overhead in all operations, because most of them 
#' require unpacking and repacking sequences, but this cost is relatively low in comparison
#' to amount of saved memory.
#' 
#' For example - \strong{dna} \strong{cln} alphabet consist of 6 values: ACGT-. They are 
#' assigned numbers 1 to 5 respectively. Those numbers in binary format take form: \code{001},
#' \code{010}, \code{011}, \code{100}, \code{101}. Each of the letters can  be coded with
#' just 3 bits instead of 8 which is demanded by \code{char} - this allows
#' us to save more than 60\% of memory spent on storage of nucleotide sequences.
#' 
#' @section tibble compatibility:
#' \code{sq} objects are compatible with \code{\link[tibble]{tibble}} class - that means one can 
#' have a \code{sq} object as a column of a \code{tibble}. There are overloaded print methods, so
#' that it is printed in pretty format.
#'
#' @name sq
NULL

#' Construct sq object from character vector
#' 
#' @description This function allows the user to construct objects of 
#' \code{\link[=sq]{class sq}} from a character vector.
#' 
#' @param sq a \code{\link{character}} vector.
#' @param type a \code{\link{character}} string indicating type of \code{sq} object that
#' is going to be constructed; supported values are "ami" for amino acid sequences,
#' "dna" for DNA sequences, "rna" for RNA sequences, "unt" for and \code{NULL} for
#' type guessing (see details)
#' @param is_clean a \code{\link{logical}} value indicating if sequences are clean.
#' or in other words - they don't contain ambiguous values; supported values are \code{TRUE} 
#' for clean sequences, \code{FALSE} for unclean sequences and \code{NULL} for auto detecting
#' (see details).
#' @param non_standard a \code{\link{character}} vector indicating non-standard letters
#' contained in sequences. If \code{NULL}, sequences will not be searched for non-standard letters
#' of length more than one. Each element of this parameter should be at least two characters 
#' long.
#' @return object of \code{\link[=sq]{class sq}} with appropriate type (one of: \strong{ami},
#' \strong{dna}, \strong{rna}, \strong{unt}, \strong{atp}).
#' 
#' @details 
#' Function \code{construct_sq} covers all possibilities of standard and non-standard types and 
#' alphabets. You can check what 'type' and 'alphabet' exactly are in \code{\link[=sq]{sq class}} 
#' documentation. Below there is a guide how function operates and how the program behaves 
#' according to the given arguments and the letters in the sequences.
#' 
#' Functions \code{construct_sq_ami}, \code{construct_sq_dna} and \code{construct_sq_rna} are
#' wrappers around \code{construct_sq} with specified \code{type} parameter - accordingly
#' "ami", "dna" or "rna". You can also pass "is_clean" parameter to those functions, but you
#' cannot pass "non_standard".
#' 
#' \code{sq} parameter should be a character vector. Each element of this vector is a biological 
#' sequence. If this parameter has length 0, object of class \code{sq} with 0 sequences will be 
#' created (if not specified, it will have \strong{dna} \strong{cln} type, which is a result of 
#' rules written below). If it contains sequences of length 0, \code{\link[=sq]{NULL}} sequences
#' will be introduced (see \emph{NULL (empty) sequences} section in \code{\link[=sq]{sq class}}).
#' 
#' \strong{Important note:} in all below cases word 'letter' stands for an element of an alphabet.
#' Letter might consist of more than one character, for example "Ala" might be a single letter.
#' However, if you want to construct or read sequences with multi-character letters, one has 
#' to specify \code{non_standard} parameter. Details of letters, alphabet and types can be 
#' found in \code{\link[=sq]{sq class}} documentation.
#' 
#' @section Simple guide to construct :
#' In most cases, just the \code{sq} parameter needs to be specified - type of sequences
#' will be guessed accordingly to rules described below. You need to pay attention, however, 
#' because for short sequences type may be guessed incorrectly - in this case you should
#' specify \code{type} and/or \code{is_clean}.
#' 
#' If your sequences contain non-standard letters, where each non-standard letter is one
#' character long, you also don't need to specify any parameter. Optionally, you can explicitly
#' do it by setting \code{type} to "unt".
#' 
#' If you want to construct sequences with multicharacter letters, you have to specify 
#' \code{non_standard} parameter, where you have to provide all non-standard letters longer
#' than one character.
#' 
#' In \code{\link[=fast-mode]{fast mode}} you have to specify both \code{type} and \code{is_clean} 
#' parameters. You cannot specify \code{non_standard} parameter in this mode. All letters
#' outside specified alphabet will be red as \code{NA} values.
#' 
#' @section Detailed guide to construct :
#' Below all possibilities that can occur during the construction of a \code{sq} object are listed.
#' 
#' In normal mode (no \code{\link[=fast-mode]{fast mode}}):
#' \itemize{
#' \item If you don't specify any other parameter than \code{sq}, function will try to guess
#' sequence type and if it is clean (it will check in exactly this order):
#' \enumerate{
#' \item If it contains only ACGT- letters, either lowercase or uppercase, type will be set
#' to \strong{dna} \strong{cln}.
#' \item If it contains only ACGU- letters, either lowercase or uppercase, type will be set
#' to \strong{dna} \strong{cln}.
#' \item If it contains any letters from 1. and 2. and additionally letters DEFHIKLMNPQRSVWY*,
#' either lowercase or uppercase, type will be set to \strong{ami} \strong{cln}.
#' \item If it contains any letters from 1. and additionally letters WSMKRYBDHVN, either
#' lowercase or uppercase, and does not contain any other letter, type will be set to 
#' \strong{dna} without \strong{cln} subtype.
#' \item If it contains any letters from 2. and additionally letters WSMKRYBDHVN, either
#' lowercase or uppercase, and does not contain any other letter, type will be set to 
#' \strong{rna} without \strong{cln} subtype.
#' \item If it contains any letters from previous points and additionally letters JOUXZ, type will
#' be set to \strong{ami} without \strong{cln} subtype.
#' \item If it contains any letters that exceed all groups mentioned above, type will be set
#' to "unt".
#' }
#' \item If you specify \code{type} parameter as "ami", "dna" or "rna" (and do not specify neither 
#' \code{is_clean} nor \code{non_standard}) type will be checked during construction: 
#' \enumerate{
#' \item If all letters in sequences fit the clean alphabet of given type, type will be set to 
#' given type with \strong{cln} subtype. 
#' \item If all letters in sequences fit the unclean alphabet of given type, type will be set to 
#' given type without \strong{cln} subtype.
#' \item If at least one sequence contains at least one letter that is not an element of an  
#' unclean alphabet of the provided type, an error will be thrown. 
#' }
#' \item If you specify both \code{type} and \code{is_clean}, function checks if letters
#' in sequences matches exactly specified alphabet (with capitalisation accuracy). If they do, 
#' type will be set to it. Otherwise, an error will be thrown.
#' \item If you specify \code{type} as "unt" and neither \code{is_clean} nor \code{non_standard},
#' type will be set to \strong{unt}. Letters won't be converted to uppercase, alphabet will
#' consist of all letters found in sequences. 
#' \item If you do not specify neither \code{type} nor \code{is_clean} and specify 
#' \code{non_standard} parameter, which should be character vector where each element is at least 
#' two characters long, all strings as specified will be detected in sequences and treated as 
#' letters in constructed \strong{atp} \code{sq}.
#' \item All other combinations of parameters are incorrect.
#' }
#' 
#' In \code{\link[=fast-mode]{fast mode}} you have to specify \code{type} (it has to have one of
#' "ami", "dna" or "rna" values) and \code{is_clean} (\code{TRUE} or \code{FALSE}). You cannot
#' specify \code{non_standard}. All letters that aren't elements of destination alphabet (with
#' a letter size accuracy) will be treated as \code{NA} values.
#' 
#' @section Handling atp and unt sequences and NA values:
#' You can convert letters into another using \code{\link{substitute_letters}} and then you
#' can use \code{\link{typify}} function to set of \code{sq} to \strong{ami}, \strong{dna}
#' or \strong{rna}. If your sequences contain \code{NA} values, use \code{\link{remove_na}}
#' 
#' @examples 
#' # saving option:
#' previous_option <- getOption("tidysq_g_fast_mode")
#' 
#' #### constructing sq in normal mode:
#' ## setting an option:
#' options(tidysq_g_fast_mode = FALSE)
#' 
#' ## constructing sq without specifying type
#' # dna cln sq
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"))
#' 
#' # rna cln sq
#' construct_sq(c("CUUAC", "UACCGGC", "GCA-ACGU"))
#' 
#' # ami cln sq
#' construct_sq(c("YQQPAVVM", "PQCFL"))
#' 
#' # ami cln sq can contain * - letter meaning end of translation:
#' construct_sq(c("MMDF*", "SYIHR*", "MGG*"))
#' 
#' # dna sq
#' construct_sq(c("TMVCCDA", "BASDT-CNN"))
#' 
#' # rna sq
#' construct_sq(c("WHDHKYN", "GCYVCYU"))
#' 
#' # ami sq
#' construct_sq(c("XYOQWWKCNJLO"))
#' 
#' # unt sq - let's assume that one wants to mark some special element in sequence with %
#' construct_sq(c("%%YAPLAA", "PLAA"))
#' 
#' ## constructing sq with type 
#' # all above examples will result in an identical if specified type as guessed
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "dna")
#' construct_sq(c("CUUAC", "UACCGGC", "GCA-ACGU"), "rna")
#' construct_sq(c("YQQPAVVM", "PQCFL"), "ami")
#' construct_sq(c("MMDF*", "SYIHR*", "MGG*"), "ami")
#' construct_sq(c("TMVCCDA", "BASDT-CNN"), "dna")
#' construct_sq(c("WHDHKYN", "GCYVCYU"), "rna")
#' construct_sq(c("XYOQWWKCNJLO"), "ami")
#' construct_sq(c("%%YAPLAA", "PLAA"), "unt")
#' 
#' # you can also use wrappers instead of parameters
#' construct_sq_dna(c("ATGC", "TCGTTA", "TT--AG"))
#' construct_sq_rna(c("CUUAC", "UACCGGC", "GCA-ACGU"))
#' construct_sq_ami(c("YQQPAVVM", "PQCFL"))
#' construct_sq_ami(c("MMDF*", "SYIHR*", "MGG*"))
#' construct_sq_dna(c("TMVCCDA", "BASDT-CNN"))
#' construct_sq_rna(c("WHDHKYN", "GCYVCYU"))
#' construct_sq_ami(c("XYOQWWKCNJLO"))
#' 
#' # One can force type other than guessed (if letters fit in the destination alphabet)
#' construct_sq(c("ATGC", "TCGTTA", "TT--AG"), "dna", is_clean = FALSE)
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
#' options(tidysq_g_fast_mode = TRUE)
#' 
#' # you cannot construct without specifying type
#' \dontrun{
#' construct_sq("CTGA")
#' }
#' construct_sq("CTGA", "dna", TRUE)
#' 
#' # you cannot construct with specifying non_standard
#' \dontrun{
#' construct_sq("mAPQ", non_standard = "mA")
#' }
#' 
#' # letters other than in specified alphabet will be treated as NA:
#' construct_sq("NCTGCNA", "dna", TRUE)
#' 
#' #### Other examples:
#' ## setting fast mode off again:
#' options(tidysq_g_fast_mode = FALSE)
#' 
#' # lowercase letters are converted to uppercase if detected type is
#' # "ami", "dna" or "rna"
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
#' options(tidysq_g_fast_mode = previous_option)
#' 
#' @seealso \code{\link{sq}} \code{\link{read_fasta}} \code{\link{tidysq-options}} 
#' \code{\link{fast-mode}} \code{\link{substitute_letters}} \code{\link{remove_na}}
#' @export
sq <- function(x,
               alphabet = guess_sq_type(x),
               NA_letter = getOption("tidysq_NA_letter"),
               safe_mode = getOption("tidysq_safe_mode")) {
  assert_character(x, any.missing = FALSE)
  assert_flag(safe_mode)
  assert_string(NA_letter)
  assert_character(alphabet, any.missing = FALSE, min.len = 1, unique = TRUE)
  
  if (length(alphabet) == 1) {
    type <- interpret_type(alphabet)
    # we can suppose that case is not important (for now)
    alphabet <- get_standard_alphabet(type)
  } else {
    type <- "atp"
  }
  alphabet <- sq_alphabet(alphabet, type)
  
  pack(x, alphabet, NA_letter, safe_mode)
}

sq_ptype <- function(str_alphabet, type)
  new_list_of(ptype = "raw",
              alphabet = sq_alphabet(str_alphabet, type),
              class = c(type_as_class(type), "sq"))
