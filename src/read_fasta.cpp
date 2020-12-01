#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::DataFrame CPP_read_fasta(const Rcpp::StringVector &file_name,
                               const Rcpp::StringVector &alphabet,
                               const Rcpp::StringVector &NA_letter,
                               const Rcpp::LogicalVector &ignore_case) {
  return export_to_R(
    read_fasta<RCPP_IT>(util::convert_to_scalar(file_name), 
                        import_alphabet_from_R(alphabet, NA_letter, ignore_case)));
}

//[[Rcpp::export]]
Rcpp::StringVector CPP_sample_fasta(const Rcpp::StringVector &file_name,
                                    const Rcpp::NumericVector &sample_size,
                                    const Rcpp::StringVector &NA_letter,
                                    const Rcpp::LogicalVector &ignore_case) {
  return export_to_R(
    sample_fasta(util::convert_to_scalar(file_name),
                 util::convert_sample_size(sample_size),
                 util::convert_to_scalar(NA_letter),
                 util::convert_to_scalar(ignore_case)));
}