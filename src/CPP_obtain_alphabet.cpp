#include "tidysq/util/obtain_alphabet.h"
#include "tidysq/tidysq-includes.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_obtain_alphabet(const Rcpp::StringVector &x,
                                       const Rcpp::IntegerVector &sample_size,
                                       const Rcpp::StringVector &NA_letter) {
    LenSq sample_size_converted = Rcpp::traits::is_infinite<INTSXP>(sample_size[0]) ? ULLONG_MAX : sample_size[0];

   return export_to_R(
           util::obtain_alphabet<RCPP_IT>(
             x, sample_size_converted, util::get_scalar_string_value(NA_letter))
           );
}

