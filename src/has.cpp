#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/ops/has.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List& x,
                            const Rcpp::StringVector& motifs) {
    return has<RCPP>(importFromR(x, "!"), Rcpp::as<std::vector<std::string>>(motifs));
}
