#include <Rcpp.h>

#include "tidysq/tidysq-includes.h"
#include "tidysq/remove_on_condition.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_remove_NA(const Rcpp::List& x,
                         const Rcpp::LogicalVector &by_letter,
                         const Rcpp::StringVector& NA_letter) {
    return export_to_R(remove_NA<RCPP_IT>(import_from_R(x, NA_letter), util::convert_to_scalar(by_letter)));
}
