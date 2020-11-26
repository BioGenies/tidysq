#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(const Rcpp::List& x,
                    const Rcpp::IntegerVector& indices,
                    const Rcpp::StringVector &NA_letter) {
    if (Rcpp::is_true(Rcpp::all(indices > 0))) {
        Rcpp::IntegerVector cpp_indices = indices - 1;
        std::pair<std::string, Sq<RCPP_IT>> ret =
                bite<RCPP_IT>(import_from_R(x, NA_letter), Rcpp::as<std::vector<long long int>>(cpp_indices));
        return Rcpp::List::create(Rcpp::Named("warning", std::get<0>(ret)),
                                  Rcpp::Named("sq", export_to_R(std::get<1>(ret))));
    } else if (Rcpp::is_true(Rcpp::all(indices < 0))) {
        Rcpp::IntegerVector cpp_indices = -indices - 1;
        Sq<RCPP_IT> ret =
                skip<RCPP_IT>(import_from_R(x, NA_letter), Rcpp::as<std::vector<long long int>>(cpp_indices));
        return Rcpp::List::create(Rcpp::Named("warning", ""),
                                  Rcpp::Named("sq", export_to_R(ret)));
    } else {
        throw std::invalid_argument("indices must be either all positive or all negative");
    }
}
