#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/ops/has.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List& x,
                            const Rcpp::StringVector& motifs,
                            const Rcpp::StringVector &NA_letter) {
    return has<RCPP>(import_from_R(x, NA_letter), 
                     Rcpp::as<std::vector<std::string>>(motifs));
}
