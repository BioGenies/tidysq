#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_remove_ambiguous(const Rcpp::List& x,
                                const Rcpp::LogicalVector &by_letter,
                                const Rcpp::StringVector& NA_letter) {
    return export_to_R(remove_ambiguous<RCPP_IT>(import_from_R(x, NA_letter), util::convert_to_scalar(by_letter)));
}
