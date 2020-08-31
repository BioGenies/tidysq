#' Set options of package
#' 
#' You can get value of an option by calling \code{getOptions(option_name)} and set its value
#' by calling \code{options(option_name = value)}, where \code{option_name} is an option name 
#' (full list of this package included below) and \code{value} is a value to assign to an option.
#' 
#' @details 
#' You can change default behavior of package using one of following options:
#' \itemize{
#' \item tidysq_a_bite_na (default "warning") - a \code{\link{character}} string specifying
#' in which way to inform user about biting sequences out of the range when using 
#' \code{\link{bite}}; possible values: "error", "warning", "message", "none",
#' \item tidysq_a_cln_sub_letters (default "warning") - a \code{\link{character}} string
#' specifying in which way to inform user
#' about dropping \code{cln} subtype of \code{sq} while using \code{\link{substitute_letters}}; 
#' possible values: "error", "warning", "message", "none",
#' \item tidysq_a_typify_small_cap_let (default "warning") - a \code{\link{character}} string
#' specifying in which way to inform user
#' about merging lowercase and uppercase letters when typifying \code{sq} object; possible 
#' values: "error", "warning", "message", "none",
#' \item tidysq_a_no_given_enc (default "warning") - a \code{\link{character}} string 
#' specifying in which way to inform user
#' about encoding unspecified letters as \code{\link{NA}} if they do appear in sequences in
#' \code{\link{encode}}; possible values: "error", "warning", "message", "none",
#' \item tidysq_g_fast_mode (default \code{FALSE}) - a \code{\link{logical}} value if to work in
#' fast mode (see 
#' \code{\link{fast-mode}})
#' \item tidysq_p_max_pillar_width (default 15) - an \code{\link{integer}} value specifying
#' pillar_shaft_sq width
#' \item tidysq_p_max_sequences (default 10) - an \code{\link{integer}} value specifying
#' maximum number of printed sequences
#' in \code{\link[=sq-print]{print sq}},
#' (see \code{\link[=sq-print]{print sq}})
#' \item tidysq_p_na_letter (default = "!") -  a \code{\link{character}} string which character
#' string to pring when 
#' \code{\link{NA}} values appear in sequences,
#' \item tidysq_p_use_color (default \code{TRUE}) - a \code{\link{logical}} value if to
#' use colorful printing
#' }
#' 
#' If value is not appropriate for given option, default value will be used.
#' 
#' @seealso \code{\link[base:option]{options}}
#' @name tidysq-options
NULL

#' Make your operations faster
#' 
#' \emph{Fast mode} is meant to improve performance of operations in package. However, turning it 
#' on is associated with less control - some of parameters checks are skipped, so user has to be
#' certain that they are correct.
#' 
#' Not all functions support \emph{fast mode} yet. Those that do not always operate in normal mode.
#' 
#' @section Turning fast mode on and off:
#' 
#' \code{Fast mode} can be turned on and off by setting package option "tidysq_g_fast_mode" to
#' accordingly \code{TRUE} and \code{FALSE} with the latter value as default. To learn more about
#' options, see \code{\link{tidysq-options}}.
#'
#' @section Functions supporting fast mode:
#' 
#' Below there is a list of functions that do support \emph{fast mode}:
#' \itemize{
#' \item \code{\link{construct_sq}}
#' \item \code{\link{read_fasta}}
#' }
#' 
#' @name fast-mode
#' @seealso \code{\link{sq}} \code{\link{tidysq-options}}
NULL