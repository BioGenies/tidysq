#' @section Motif capabilities and restrictions:
#' There are more options than to simply create a motif that is a string
#' representation of searched subsequence. For example, when using this function
#' with any of standard types, i.e. \strong{ami}, \strong{dna} or \strong{rna},
#' the user can create a motif with ambiguous letters. In this case the engine
#' will try to match any of possible meanings of this letter. For example, take
#' "B" from extended DNA alphabet. It means "not A", so it can be matched with
#' "C", "G" and "T", but also "B", "Y" (either "C" or "T"), "K" (either "G" or
#' "T") and "S" (either "C" or "G").
#'
#' Full list of ambiguous letters with their meaning can be found on IUPAC site.
#'
#' Motifs are also restricted in that the alphabets of \code{sq} objects on
#' which search operations are conducted cannot contain "^" and "$" symbols.
#' These two have a special meaning - they are used to indicate beginning and
#' end of sequence respectively and can be used to limit the position of matched
#' subsequences.
