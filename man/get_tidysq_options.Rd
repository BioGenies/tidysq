% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util-options.R
\name{get_tidysq_options}
\alias{get_tidysq_options}
\alias{tidysq-options}
\title{Obtain current state of tidysq options}
\usage{
get_tidysq_options()
}
\value{
A \code{\link[=namedList-class]{named list}} with selected option
values.
}
\description{
Subsets all global options to display those related to
\pkg{tidysq} package.
}
\details{
The user can display value of selected option by calling
\code{getOptions(option_name)} and set its value with
\code{options(option_name = value)}, where \code{option_name} is an option
name and \code{value} is a value to assign to an option.

Full list of options included in \pkg{tidysq} package is listed below:
\itemize{
\item tidysq_NA_letter [\code{character(1)}]\cr
 A letter to be used when printing, constructing or interpreting \code{NA}
 value. Defaults to \code{"!"}.
\item tidysq_on_warning [\code{"silent" || "message" || "warning" || "error"}]\cr
 Determines the method of handling warning message. Setting \code{"error"}
 makes any warning throw an exception and stop execution of the code. The
 difference between \code{"message"} and \code{"warning"} is that while both
 display warning text to the console, only the latter registers it so that
 it can be accessed with a call to \code{warnings()}. Lastly, \code{"silent"}
 setting causes any warnings to be completely ignored. Default value is
 \code{"warning"}.
\item tidysq_pillar_max_width [\code{integer(1)}]\cr
 Determines max width of a column of \code{sq} class within a
 \code{\link[tibble]{tibble}}. Default value is 15.
\item tidysq_print_max_sequences [\code{integer(1)}]\cr
 Controls maximum number of sequences printed to console. If an \code{sq}
 object is longer than this value, then only first
 \code{tidysq_print_max_sequences} are printed, just like in any R vector.
 Default value is 10.
\item tidysq_print_use_color [\code{logical(1)}]\cr
 Determines whether coloring should be used to increase readability of text
 printed to console. While it is advised to keep this option turned on due
 to above concern, some environments may not support coloring and thus
 turning it off can be necessary. Defaults to \code{TRUE}.
\item tidysq_safe_mode [\code{logical(1)}]\cr
 Default value is \code{FALSE}. When turned on, safe mode guarantees that
 \code{NA} appears within a sequence if and only if input sequence contains
 value passed with \code{NA_letter}. This means that resulting type might be
 different to the one passed as argument, if there are letters in a sequence
 that does not appear in the original alphabet.
}
}
\seealso{
Functions that display sequence info:
\code{\link{sqprint}}
}
\concept{display_functions}
