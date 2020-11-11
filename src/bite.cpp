#include <Rcpp.h>

#include "tidysq/tidysq-includes.h"
#include "tidysq/bite.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(const Rcpp::List& x,
                    const Rcpp::IntegerVector& indices,
                    const Rcpp::StringVector &NA_letter) {
    auto ret = bite<RCPP>(import_from_R(x, NA_letter), Rcpp::as<std::vector<int>>(indices));
    return Rcpp::List::create(Rcpp::Named("warning", std::get<0>(ret)),
                              Rcpp::Named("sq", export_to_R(std::get<1>(ret))));
}
