#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_obtain_alphabet(const Rcpp::StringVector &x,
                                       const Rcpp::NumericVector &sample_size,
                                       const tidysq::Letter &NA_letter,
                                       const bool &ignore_case) {
    return export_to_R(
           util::obtain_alphabet<RCPP_IT>(x,
                                          util::convert_sample_size(sample_size),
                                          NA_letter,
                                          ignore_case)
           );
}

