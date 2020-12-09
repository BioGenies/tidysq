#include "tidysq/Rcpp-import.h"
#include "tidysq/ops/apply_R_function.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_apply_R_function(const Rcpp::List &x,
                                const Rcpp::Function &fun,
                                const Rcpp::LogicalVector &single_string,
                                const Rcpp::StringVector &NA_letter) {
    const bool single_string_converted = util::convert_to_scalar(single_string);
    if (single_string_converted) {
        return apply_R_function<RCPP_IT, STRING_PT>(import_from_R(x, NA_letter), fun);
    } else {
        return apply_R_function<RCPP_IT, STRINGS_PT>(import_from_R(x, NA_letter), fun);
    }
}