#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/ops/bite.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(const Rcpp::List& x, const Rcpp::IntegerVector& indices) {
    auto ret = bite<RCPP>(importFromR(x, "!"), Rcpp::as<std::vector<int>>(indices));
    return Rcpp::List::create(Rcpp::Named("warning") = std::get<0>(ret),
                              Rcpp::Named("sq") = std::get<1>(ret).exportToR());
}
