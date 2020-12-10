#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::DataFrame CPP_read_fasta(const std::string &file_name,
                               const Rcpp::StringVector &alphabet,
                               const tidysq::Letter &NA_letter,
                               const bool &ignore_case) {
  return export_to_R(
    io::read_fasta<RCPP_IT>(file_name,
                            import_alphabet_from_R(alphabet, NA_letter, ignore_case)));
}

//[[Rcpp::export]]
Rcpp::StringVector CPP_sample_fasta(const std::string &file_name,
                                    const Rcpp::NumericVector &sample_size,
                                    const tidysq::Letter &NA_letter,
                                    const bool &ignore_case) {
  return export_to_R(
    io::sample_fasta(file_name,
                     util::convert_sample_size(sample_size),
                     NA_letter,
                     ignore_case));
}