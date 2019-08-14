#' Set options of package
#' 
#' @usage 
#' \code{getOptions(option_name)}
#' 
#' \code{options(option_name = value)}
#' 
#' @param option_name an option name to set or get
#' @param value a value to assing to an option
#' 
#' @return in case of \code{getOptions} - value of the option
#' 
#' @details 
#' You can change default behaviour of package using one of following options:
#' \itemize{
#' \item tidysq_bite_na_action (default "warning") - a \code{\link{character}} string specifying 
#' in which way to inform user about biting sequences out of the range when using 
#' \code{\link{bite}}; possible values: "error", "warning", "message", "none",
#' \item tidysq_subsitute_letters_cln (default "warning") - a \code{\link{character}} string 
#' specifying in which way to inform user
#' about droping \code{cln} subtype of \code{sq} while using \code{\link{substitute_letters}}; 
#' possible values: "error", "warning", "message", "none",
#' \item tidysq_typify_small_cap_let (default "warning") - a \code{\link{character}} string 
#' specifying in which way to inform user
#' about merging lowercase and uppercase letters when typifying \code{sq} object; possible 
#' values: "error", "warning", "message", "none",
#' \item tidysq_encode_no_given_action (default "warning") - a \code{\link{character}} string 
#' specifying in which way to inform user
#' about encoding unspecified letters as \code{\link{NA}} if they do appear in sequences in
#' \code{\link{encode}}; possible values: "error", "warning", "message", "none",
#' \item tidysq_max_pillar_sq_width (default 15) - an \code{\link{integer}} value specyfying 
#' pillar_shaft_sq width
#' \item tidysq_max_print_sequences (default 10) - an \code{\link{integer}} value specyfying 
#' maximum number of printed sequences
#' in \code{\link[sq:print.sq]{print sq}},
#' \item tidysq_colorful_sq_print (default \code{TRUE}) - a \code{\link{logical}} value if to 
#' use colorful printing
#' (see \code{\link[sq:print.sq]{print sq}})
#' \item tidysq_na_print_char (default = "!") -  a \code{\link{character}} string which character 
#' string to pring when 
#' \code{\link{NA}} values appear in sequences,
#' \item tidysq_fast_mode (default \code{FALSE}) - a \code{\link{logical}} value if to work in 
#' fast mode (see 
#' \code{\link{fast-mode}})
#' }
#' 
#' If value is not appropriate for given option, default value will be used.
#' 
#' @seealso \code{\link[base:option]{options}}
#' @name sq-options
NULL