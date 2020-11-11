#include "tidysq/util/obtain_alphabet.h"
#include "tidysq/tidysq-includes.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_obtain_alphabet(const Rcpp::StringVector &x,
                                       const Rcpp::IntegerVector &sample_size,
                                       const Rcpp::StringVector &NA_letter) {
   return export_to_R(
           util::obtain_alphabet<RCPP_IT>(
             x, sample_size[0], util::get_scalar_string_value(NA_letter))
           );
}

