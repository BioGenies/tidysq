% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_fasta.R
\name{write_fasta}
\alias{write_fasta}
\alias{write_fasta.sq}
\alias{write_fasta.data.frame}
\title{Save sq to fasta file}
\usage{
write_fasta(x, ...)

\method{write_fasta}{sq}(
  x,
  name,
  file,
  width = 80,
  NA_letter = getOption("tidysq_NA_letter"),
  ...
)

\method{write_fasta}{data.frame}(
  x,
  file,
  .sq = "sq",
  .name = "name",
  width = 80,
  NA_letter = getOption("tidysq_NA_letter"),
  ...
)
}
\arguments{
\item{x}{[\code{sq}]\cr
An object this function is applied to.}

\item{...}{further arguments to be passed from or to other methods.}

\item{name}{[\code{character}]\cr
Vector of sequence names. Must be of the same length as \code{sq} object.}

\item{file}{[\code{character(1)}]\cr
Absolute path to file to write to.}

\item{width}{[\code{integer(1)}]\cr
Maximum number of characters to put in each line of file. Must be positive.}

\item{NA_letter}{[\code{character(1)}]\cr
A string that is used to interpret and display \code{NA} value in the
context of \code{\link[=sq-class]{sq class}}. Default value equals to
"\code{!}".}

\item{.sq}{[\code{character(1)}]\cr
Name of a column that stores sequences.}

\item{.name}{[\code{character(1)}]\cr
Name of a column that stores names (unique identifiers).}
}
\value{
No value is returned.
}
\description{
Writes \code{\link[=sq-class]{sq}} objects with their names to
a fasta file.
}
\details{
Whenever a name has more letters than \code{width} parameter, nothing
happens, as only sequences are split to fit within designated space.
}
\examples{
\dontrun{
sq_dna <- sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"),
             alphabet = "dna_bsc")
write_fasta(sq_dna,
            c("bat", "cat", "rat", "elephant_swallowed_by_A_snake"),
            tempfile())
}

# It can be a part of tidyverse pipeline:
library(dplyr)
fasta_file <- system.file(package = "tidysq", "examples/example_aa.fasta")
read_fasta(fasta_file) \%>\%
  mutate(name = toupper(name)) \%>\%
  write_fasta(tempfile())

}
\seealso{
Functions from output module:
\code{\link{as.character.sq}()},
\code{\link{as.matrix.sq}()},
\code{\link{as.sq}()},
\code{\link{export_sq}()}
}
\concept{output_functions}
