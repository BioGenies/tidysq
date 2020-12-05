#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_reverse(const Rcpp::List& x,
                       const tidysq::Letter& NA_letter) {
    return export_to_R(reverse<RCPP_IT>(import_from_R(x, NA_letter)));
}
