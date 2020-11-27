#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_obtain_alphabet(const Rcpp::StringVector &x,
                                       const Rcpp::NumericVector &sample_size,
                                       const Rcpp::StringVector &NA_letter,
                                       const Rcpp::LogicalVector &ignore_case) {
    LenSq sample_size_converted = Rcpp::traits::is_infinite<REALSXP>(sample_size[0]) ? ULLONG_MAX : sample_size[0];

   return export_to_R(
           util::obtain_alphabet<RCPP_IT>(
             x, sample_size_converted, util::convert_to_scalar(NA_letter), util::convert_to_scalar(ignore_case))
           );
}

