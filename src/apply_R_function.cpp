#include "tidysq/Rcpp-import.h"
#include "tidysq/ops/apply_R_function.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_apply_R_function(const Rcpp::List &x,
                                const Rcpp::Function &fun,
                                const bool &single_string,
                                const tidysq::Letter &NA_letter) {
    if (single_string) {
        return apply_R_function<RCPP_IT, STRING_PT>(import_from_R(x, NA_letter), fun);
    } else {
        return apply_R_function<RCPP_IT, STRINGS_PT>(import_from_R(x, NA_letter), fun);
    }
}