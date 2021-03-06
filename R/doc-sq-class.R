#' sq: class for keeping biological sequences tidy
#' 
#' An object of class \strong{sq} represents a list of biological sequences. It
#' is the main internal format of the \pkg{tidysq} package and most functions
#' operate on it. The storage method is memory-optimized so that objects
#' require as little memory as possible (details below).
#' 
#' @section Construction/reading/import of sq objects:
#' There are multiple ways of obtaining \code{sq} objects:
#' \itemize{
#' \item constructing from a \code{\link{character}} vector with
#' \code{\link{sq}} method,
#' \item constructing from another object with \code{\link{as.sq}} method,
#' \item reading from the FASTA file with \code{\link{read_fasta}},
#' \item importing from a format of other package like \pkg{ape} or
#' \pkg{Biostrings} with \code{\link{import_sq}}.
#' }
#'
#' \strong{Important note:} A manual assignment of a class \code{sq} to an
#' object is \strong{strongly discouraged} - due to the usage of low-level
#' functions for bit packing such assignment may lead to calling one of those
#' functions during operating on object or even printing it which can cause
#' a crash of R session and, in consequence, loss of data.
#'
#' @section Export/writing of sq objects:
#' There are multiple ways of saving \code{sq} objects or converting them into
#' other formats:
#' \itemize{
#' \item converting into a character vector with
#' \code{\link[=as.character.sq]{as.character}} method,
#' \item converting into a character matrix with
#' \code{\link[=as.matrix.sq]{as.matrix}} method,
#' \item saving as FASTA file with \code{\link{write_fasta}},
#' \item exporting into a format of other package like \code{ape} or
#' \code{Biostrings} with \code{\link{export_sq}}.
#' }
#'
#' @section Ambiguous letters:
#' This package is meant to handle amino acid, DNA and RNA sequences. IUPAC
#' standard for one letter codes includes ambiguous bases that are used to
#' describe more than one basic standard base. For example, "\code{B}" in the
#' context of DNA code means "any of C, G or T". As there are operations that
#' make sense only for unambiguous bases (like \code{\link{translate}}), this
#' package has separate types for sequences with "basic" and "extended"
#' alphabet.
#'
#' @section Types of sq:
#' There is need to differentiate \code{sq} objects that keep different types
#' of sequences (DNA, RNA, amino acid), as they use different alphabets.
#' Furthermore, there are special types for handling non-standard sequence
#' formats.
#' 
#' Each \strong{sq} object has exactly one of \strong{types}:
#' \itemize{
#' \item \strong{ami_bsc} - (\emph{amino acids}) represents a list of sequences
#' of amino acids (peptides or proteins),
#' \item \strong{ami_ext} - same as above, but with possible usage of
#' ambiguous letters,
#' \item \strong{dna_bsc} - (\emph{DNA}) represents a list of DNA sequences,
#' \item \strong{dna_ext} - same as above, but with possible usage of
#' ambiguous letters,
#' \item \strong{rna_bsc} - (\emph{RNA}) represents a list of RNA sequences
#' (together with DNA above often collectively called "nucleotide sequences"),
#' \item \strong{rna_ext} - same as above, but with possible usage of
#' ambiguous letters,
#' \item \strong{unt} - (\emph{untyped}) represents a list of sequences that do
#' not have specified type. They are mainly result of reading sequences from
#' a file that contains some letters that are not in standard nucleotide or
#' amino acid alphabets and user has not specified them explicitly. They should
#' be converted to other \strong{sq} classes (using functions like
#' \code{\link{substitute_letters}} or \code{\link{typify}}),
#' \item \strong{atp} - (\emph{atypical}) represents sequences that have an
#' alphabet different from standard alphabets - similarly to \strong{unt}, but
#' user has been explicitly informed about it. They are result of constructing
#' sequences or reading from file with provided custom alphabet (for details
#' see \code{\link{read_fasta}} and \code{\link{sq}} function). They are also
#' result of using function \code{\link{substitute_letters}} - users can use
#' it to for example simplify an alphabet and replace several letters by one.
#' }
#'
#' For clarity, \strong{ami_bsc} and \strong{ami_ext} types are often referred
#' to collectively as \strong{ami} when there is no need to explicitly specify
#' every possible type. The same applies to \strong{dna} and \strong{rna}.
#'
#' \code{sq} object type is printed when using overloaded method
#' \code{\link[=sq-print]{print}}. It can be also checked and obtained as
#' a value (that may be passed as argument to function) by using
#' \code{\link{sq_type}}.
#'
#' @section Alphabet:
#' See \code{\link{alphabet}}.
#'
#' The user can obtain an alphabet of the \code{sq} object using the
#' \code{\link{alphabet}} function. The user can check which letters are
#' invalid (i.e. not represented in standard amino acid or nucleotide
#' alphabet) in each sequence of given \code{sq} object by using
#' \code{\link{find_invalid_letters}}. To substitute one letter with another
#' use \code{\link{substitute_letters}}.
#'
#' @section Missing/Not Available values:
#' There is a possibility of introducing \code{\link{NA}} values into
#' sequences. \code{NA} value does not represents gap (which are represented by
#' "\code{-}") or wildcard elements ("\code{N}" in the case of nucleotides and
#' "\code{X}" in the case of amino acids), but is used as a representation of
#' an empty position or invalid letters (not represented in nucleotide or amino
#' acid alphabet).
#'
#' \code{NA} does not belong to any alphabet. It is printed as "\code{!}" and,
#' thus, it is highly unrecommended to use "\code{!}" as special letter in
#' \strong{atp} sequences (but print character can be changed in options, see
#' \code{\link{tidysq-options}}).
#'
#' \code{NA} might be introduced by:
#' \itemize{
#' \item reading fasta file with non-standard letters with
#' \code{\link{read_fasta}} with \code{safe_mode} argument set to \code{TRUE},
#' \item replacing a letter with \code{NA} value with
#' \code{\link{substitute_letters}},
#' \item subsetting sequences beyond their lengths with \code{\link{bite}}.
#' }
#'
#' The user can convert sequences that contain \code{NA} values into
#' \code{NULL} sequences with \code{\link{remove_na}}.
#' 
#' @section NULL (empty) sequences:
#' \code{NULL} sequence is a sequence of length 0.
#'
#' \code{NULL} sequences might be introduced by:
#' \itemize{
#' \item constructing \code{sq} object from character string of length zero,
#' \item using the \code{\link{remove_ambiguous}} function,
#' \item using the \code{\link{remove_na}} function,
#' \item subsetting \code{sq} object with \code{\link{bite}} function (and
#' negative indices that span at least \code{-1:-length(sequence)}.
#' }
#'
#' @section Storage format:
#' \code{sq} object is, in fact, \strong{list of raw vectors}. The fact that it
#' is list implies that the user can concatenate \code{sq} objects using
#' \code{\link[=sq-concatenate]{c}} method and subset them using
#' \code{\link[=sq-extract]{extract operator}}. Alphabet is kept as an
#' attribute of the object.
#' 
#' Raw vectors are the most efficient way of storage - each letter of a
#' sequence is assigned an integer (its index in alphabet of \code{sq} object).
#' Those integers in binary format fit in less than 8 bits, but normally are
#' stored on 16 bits. However, thanks to bit packing it is possible to remove
#' unused bits and store numbers more tightly. This means that all operations
#' must either be implemented with this packing in mind or accept a little time
#' overhead induced by unpacking and repacking sequences. However, this cost
#' is relatively low in comparison to amount of saved memory.
#' 
#' For example - \strong{dna_bsc} alphabet consists of 5 values: ACGT-. They
#' are assigned numbers 0 to 4 respectively. Those numbers in binary format
#' take form: \code{000}, \code{001}, \code{010}, \code{011}, \code{100}. Each
#' of these letters can be coded with just 3 bits instead of 8 which is
#' demanded by \code{char} - this allows us to save more than 60\% of memory
#' spent on storage of basic nucleotide sequences.
#' 
#' @section tibble compatibility:
#' \code{sq} objects are compatible with \code{\link[tibble]{tibble}} class -
#' that means one can have an \code{sq} object as a column of a \code{tibble}.
#' There are overloaded print methods, so that it is printed in pretty format.
#'
#' @name sq-class
NULL
