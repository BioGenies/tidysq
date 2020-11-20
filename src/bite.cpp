#include <Rcpp.h>

#include "tidysq/tidysq-includes.h"
#include "tidysq/bite.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(const Rcpp::List& x,
                    const Rcpp::IntegerVector& indices,
                    const Rcpp::StringVector &NA_letter) {
    Rcpp::IntegerVector cpp_indices;
    if (Rcpp::is_true(Rcpp::all(indices > 0))) {
        cpp_indices = indices - 1;
    } else if (Rcpp::is_true(Rcpp::all(indices < 0))) {
        cpp_indices = indices + 1;
    } else {
        throw std::invalid_argument("indices must be either all positive or all negative");
    }
    auto ret = bite<RCPP_IT>(import_from_R(x, NA_letter), Rcpp::as<std::vector<int>>(cpp_indices));
    return Rcpp::List::create(Rcpp::Named("warning", std::get<0>(ret)),
                              Rcpp::Named("sq", export_to_R(std::get<1>(ret))));
}
